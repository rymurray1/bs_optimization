"""
Copper-Focused Mining Process Optimization Model
Predicts optimal process parameters based on ore composition data
"""

import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.multioutput import MultiOutputRegressor
import matplotlib.pyplot as plt
import seaborn as sns
from typing import List, Dict, Tuple
import warnings
warnings.filterwarnings('ignore')


class MiningDataLoader:
    """Load and parse composition and process data from stacked CSV format"""

    def __init__(self, composition_file: str, process_file: str):
        self.composition_file = composition_file
        self.process_file = process_file
        self.composition_datasets = []
        self.process_datasets = []

    def parse_stacked_csv(self, filepath: str) -> List[pd.DataFrame]:
        """Parse CSV file containing multiple stacked tables separated by blank lines"""
        datasets = []
        current_table = []
        current_header = None

        with open(filepath, 'r', encoding='utf-8') as f:
            lines = f.readlines()

        for line in lines:
            line = line.strip()

            # Skip empty lines - they separate datasets
            if not line:
                if current_table and current_header:
                    # Save the current dataset
                    df = pd.DataFrame(current_table, columns=current_header)
                    datasets.append(df)
                    current_table = []
                    current_header = None
                continue

            # Split by comma
            values = [v.strip() for v in line.split(',')]

            # Check if this looks like a header line
            # Headers contain "Composite" or "Deposit" or "Criteria" as first column
            # and units like "(g/t)", "(%)","Units" etc
            first_val = values[0] if values else ""
            is_header = (
                first_val in ['Composite', 'Deposit', 'Criteria', 'Parameter'] or
                any('(' in v and ')' in v for v in values[:3]) or  # Contains units
                'Units' in values
            )

            if is_header and current_header is None:
                # This is a header line
                current_header = values
            elif current_header is not None:
                # This is a data row
                current_table.append(values)
            # else: skip metadata lines (dates, project names, etc)

        # Don't forget the last dataset
        if current_table and current_header:
            df = pd.DataFrame(current_table, columns=current_header)
            datasets.append(df)

        return datasets

    def clean_numeric_value(self, value: str) -> float:
        """Convert string value to numeric, handling multi-values and special cases"""
        if pd.isna(value) or value == '' or value == '-':
            return np.nan

        # Handle less-than symbols (<5 becomes 5)
        if '<' in str(value):
            return float(value.replace('<', ''))

        # Handle multiple values (take the average)
        if ' ' in str(value) and len(str(value).split()) > 1:
            nums = []
            for v in str(value).split():
                try:
                    nums.append(float(v))
                except:
                    pass
            if nums:
                return np.mean(nums)
            return np.nan

        # Standard numeric conversion
        try:
            return float(value)
        except:
            return np.nan

    def load_data(self) -> Tuple[List[pd.DataFrame], List[pd.DataFrame]]:
        """Load both composition and process datasets"""
        print("Loading composition data...")
        self.composition_datasets = self.parse_stacked_csv(self.composition_file)
        print(f"Found {len(self.composition_datasets)} composition datasets")

        print("Loading process data...")
        self.process_datasets = self.parse_stacked_csv(self.process_file)
        print(f"Found {len(self.process_datasets)} process datasets")

        return self.composition_datasets, self.process_datasets


class FeatureEngineer:
    """Engineer features relevant to copper optimization"""

    @staticmethod
    def add_cu_features(df: pd.DataFrame) -> pd.DataFrame:
        """Add copper-specific derived features"""
        df = df.copy()

        # Identify metal columns (convert to numeric first)
        metal_cols = ['Cu', 'Au', 'Ag', 'Fe', 'S', 'Pb', 'Zn', 'Mn', 'As', 'Sb']

        # Ensure numeric
        for col in metal_cols:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')

        # Cyanide consumers (Cu, Fe, S compete with Au/Ag for cyanide)
        if all(c in df.columns for c in ['Cu', 'Fe', 'S']):
            df['CN_Consumers'] = df['Cu'] + df['Fe'] + df['S']

        # Cu intensity relative to precious metals
        if all(c in df.columns for c in ['Cu', 'Au', 'Ag']):
            df['Cu_to_PreciousMetals'] = df['Cu'] / (df['Au'] + df['Ag'] + 0.001)  # Avoid div by 0

        # Total sulfides (affects leaching)
        sulfide_cols = [c for c in ['S', 'Cu', 'Pb', 'Zn'] if c in df.columns]
        if sulfide_cols:
            df['Total_Sulfides'] = df[sulfide_cols].sum(axis=1)

        # Base metal content
        base_metal_cols = [c for c in ['Cu', 'Pb', 'Zn'] if c in df.columns]
        if base_metal_cols:
            df['Total_BaseMetals'] = df[base_metal_cols].sum(axis=1)

        # Deleterious elements (Mn, As, Sb)
        deleterious_cols = [c for c in ['Mn', 'As', 'Sb'] if c in df.columns]
        if deleterious_cols:
            df['Deleterious_Elements'] = df[deleterious_cols].sum(axis=1)

        return df


class CuOptimizationModel:
    """Machine learning model for predicting process parameters with Cu focus"""

    def __init__(self, model_type='random_forest'):
        self.model_type = model_type
        self.model = None
        self.scaler_X = StandardScaler()
        self.scaler_y = StandardScaler()
        self.feature_names = None
        self.target_names = None
        self.feature_importance = None

    def prepare_data(self, X: pd.DataFrame, y: pd.DataFrame) -> Tuple[np.ndarray, np.ndarray]:
        """Prepare and scale data for modeling"""
        # Store names
        self.feature_names = X.columns.tolist()
        self.target_names = y.columns.tolist()

        # Convert to numeric and handle missing values
        X_clean = X.apply(pd.to_numeric, errors='coerce').fillna(X.median())
        y_clean = y.apply(pd.to_numeric, errors='coerce').fillna(y.median())

        # Scale features
        X_scaled = self.scaler_X.fit_transform(X_clean)
        y_scaled = self.scaler_y.fit_transform(y_clean)

        return X_scaled, y_scaled

    def train(self, X: pd.DataFrame, y: pd.DataFrame):
        """Train the model"""
        X_scaled, y_scaled = self.prepare_data(X, y)

        # Choose base estimator
        if self.model_type == 'random_forest':
            base_model = RandomForestRegressor(
                n_estimators=100,
                max_depth=10,
                min_samples_split=3,
                random_state=42,
                n_jobs=-1
            )
        elif self.model_type == 'gradient_boosting':
            base_model = GradientBoostingRegressor(
                n_estimators=100,
                max_depth=5,
                learning_rate=0.1,
                random_state=42
            )
        else:
            raise ValueError(f"Unknown model type: {self.model_type}")

        # Use MultiOutputRegressor for multiple target variables
        self.model = MultiOutputRegressor(base_model)

        print(f"Training {self.model_type} model...")
        self.model.fit(X_scaled, y_scaled)
        print("Training complete!")

        # Calculate feature importance (for Random Forest)
        if self.model_type == 'random_forest':
            importances = []
            for estimator in self.model.estimators_:
                importances.append(estimator.feature_importances_)
            self.feature_importance = pd.DataFrame({
                'Feature': self.feature_names,
                'Importance': np.mean(importances, axis=0)
            }).sort_values('Importance', ascending=False)

    def predict(self, X: pd.DataFrame) -> pd.DataFrame:
        """Make predictions for new compositions"""
        X_clean = X[self.feature_names].apply(pd.to_numeric, errors='coerce').fillna(X.median())
        X_scaled = self.scaler_X.transform(X_clean)

        y_pred_scaled = self.model.predict(X_scaled)
        y_pred = self.scaler_y.inverse_transform(y_pred_scaled)

        return pd.DataFrame(y_pred, columns=self.target_names, index=X.index)

    def evaluate(self, X: pd.DataFrame, y: pd.DataFrame) -> Dict:
        """Evaluate model performance with cross-validation"""
        X_scaled, y_scaled = self.prepare_data(X, y)

        # Cross-validation scores
        cv_scores = cross_val_score(
            self.model,
            X_scaled,
            y_scaled,
            cv=5,
            scoring='r2',
            n_jobs=-1
        )

        results = {
            'cv_mean_r2': cv_scores.mean(),
            'cv_std_r2': cv_scores.std(),
            'cv_scores': cv_scores
        }

        print(f"\nModel Evaluation ({self.model_type}):")
        print(f"Cross-validated R² Score: {results['cv_mean_r2']:.3f} (+/- {results['cv_std_r2']:.3f})")

        return results

    def show_feature_importance(self, top_n=15):
        """Display top feature importances"""
        if self.feature_importance is not None:
            print(f"\nTop {top_n} Most Important Features:")
            print(self.feature_importance.head(top_n).to_string(index=False))
            return self.feature_importance.head(top_n)
        else:
            print("Feature importance not available for this model type")
            return None


def main():
    """Main execution function"""

    # File paths
    composition_file = r"C:\Users\Ryan Murray\Desktop\main_model\main\composition_data.csv"
    process_file = r"C:\Users\Ryan Murray\Desktop\main_model\main\process_data.csv"

    # Load data
    loader = MiningDataLoader(composition_file, process_file)
    comp_datasets, proc_datasets = loader.load_data()

    # For now, let's analyze the first dataset as an example
    print("\n" + "="*80)
    print("DATASET STRUCTURE ANALYSIS")
    print("="*80)

    if comp_datasets:
        print("\nFirst Composition Dataset:")
        print(comp_datasets[0].head())
        print(f"\nShape: {comp_datasets[0].shape}")
        print(f"Columns: {comp_datasets[0].columns.tolist()}")

    if proc_datasets:
        print("\nFirst Process Dataset:")
        print(proc_datasets[0].head(10))
        print(f"\nShape: {proc_datasets[0].shape}")
        print(f"Columns: {proc_datasets[0].columns.tolist()}")

    print("\n" + "="*80)
    print("Note: The data structure needs further processing to align composition")
    print("and process parameters. Each dataset may have different formats.")
    print("="*80)

    # TODO: Next steps
    print("\nNext Steps:")
    print("1. Standardize dataset formats")
    print("2. Align composition rows with process parameter rows")
    print("3. Extract numeric values from multi-value cells")
    print("4. Create feature matrix (X) and target matrix (y)")
    print("5. Train and validate model")


if __name__ == "__main__":
    main()
