"""
Optimization model for flotation recovery using scipy.optimize.

This model takes a given input feed composition and optimizes flotation recovery
by adjusting key process parameters such as:
- Specific power (sp_power)
- Specific gas rate (sp_gas_rate)
- Frother concentration (frother_conc)
- Retention time (ret_time)
- Number of cells (num_cells)
"""

import numpy as np
from scipy.optimize import minimize, differential_evolution
from physical_parameters import flot_rec, grind_break
import warnings
warnings.filterwarnings('ignore')


class FlotationOptimizer:
    """
    Optimizer for flotation recovery rates based on physical parameters.
    """

    def __init__(self, feed_composition, constant_params):
        """
        Initialize the optimizer with feed composition and constant parameters.

        Args:
            feed_composition: Dict containing feed mass distribution by size and grade
            constant_params: Dict of truly constant flotation parameters (not optimized)
        """
        self.feed_composition = feed_composition
        self.constant_params = constant_params

    def objective_function(self, x):
        """
        Objective function to maximize copper recovery (minimize negative copper recovery).

        Args:
            x: Array of decision variables:
               [sp_power, sp_gas_rate, frother_conc, ret_time, num_cells,
                air_fraction, slurry_fraction, bubble_z_pot, cell_volume, froth_height]

        Returns:
            float: Negative copper recovery (for minimization)
        """
        # Unpack decision variables
        sp_power, sp_gas_rate, frother_conc, ret_time, num_cells, \
        air_fraction, slurry_fraction, bubble_z_pot, cell_volume, froth_height = x

        num_cells = int(round(num_cells))  # Must be integer

        # Get feed data
        size_classes = self.feed_composition['size_classes']
        grade_classes = self.feed_composition['grade_classes']
        throughput_distribution = self.feed_composition['throughput_distribution']
        sg_by_grade = self.feed_composition['sg_by_grade']
        contact_angles_by_grade = self.feed_composition['contact_angles_by_grade']
        particle_z_pot_by_grade = self.feed_composition['particle_z_pot_by_grade']
        grade_by_grade = self.feed_composition['grade_by_grade']

        # Constant parameters (not optimized)
        frother = self.constant_params['frother']
        permitivity = self.constant_params['permitivity']
        dielectric = self.constant_params['dielectric']
        pe = self.constant_params['pe']
        water_or_particle = self.constant_params['water_or_particle']
        dbl_cell_area = self.constant_params['dbl_cell_area']
        dbl_bbl_ratio = self.constant_params['dbl_bbl_ratio']

        total_feed_cu = 0
        total_recovered_cu = 0

        try:
            # Calculate copper recovery for each size-grade combination
            for size_idx, particle_size in enumerate(size_classes):
                for grade_idx, grade_pct in enumerate(grade_classes):
                    contact_angle_deg = contact_angles_by_grade[grade_idx]
                    sg_for_grade = sg_by_grade[grade_idx]
                    particle_z_pot_for_grade = particle_z_pot_by_grade[grade_idx]
                    grade_value = grade_by_grade[grade_idx]

                    mass_in_class = throughput_distribution[size_idx][grade_idx]
                    cu_in_class = mass_in_class * grade_value / 100  # Copper mass in this class

                    # Calculate recovery
                    recovery = flot_rec(
                        sg_for_grade, sp_power, sp_gas_rate, air_fraction, slurry_fraction,
                        particle_z_pot_for_grade, bubble_z_pot, num_cells, ret_time, cell_volume,
                        froth_height, frother, frother_conc, particle_size, grade_value,
                        contact_angle_deg, permitivity, dielectric, pe, water_or_particle,
                        dbl_cell_area, dbl_bbl_ratio
                    )

                    total_feed_cu += cu_in_class
                    total_recovered_cu += cu_in_class * recovery

            overall_cu_recovery = total_recovered_cu / total_feed_cu if total_feed_cu > 0 else 0

            # Return negative copper recovery for minimization
            return -overall_cu_recovery

        except Exception as e:
            # Return large penalty if calculation fails
            return 1e6

    def optimize(self, method='differential_evolution', bounds=None):
        """
        Optimize flotation recovery.

        Args:
            method: Optimization method ('differential_evolution', 'slsqp', 'trust-constr')
            bounds: Tuple of (min, max) bounds for each variable
                   [(sp_power_min, sp_power_max), (sp_gas_rate_min, sp_gas_rate_max), ...]

        Returns:
            dict: Optimization results including optimal parameters and recovery
        """
        # Default bounds if not provided
        if bounds is None:
            bounds = [
                (0.8, 1.3),      # sp_power (kW/m^3)
                (0.8, 2.0),      # sp_gas_rate (cm/s)
                (10, 100),       # frother_conc (mg/L)
                (10.0, 30.0),    # ret_time (minutes)
                (5, 15),         # num_cells
                (0.05, 0.15),    # air_fraction
                (0.10, 0.30),    # slurry_fraction
                (-0.05, -0.01),  # bubble_z_pot (V)
                (50, 300),       # cell_volume (m^3)
                (0.10, 0.30)     # froth_height (m, converted from 10-30 cm)
            ]

        # Initial guess (current parameters)
        x0 = [1.0, 1.5, 50, 20.0, 10, 0.10, 0.20, -0.03, 175, 0.165]

        print(f"\nStarting optimization using {method}...")
        print(f"Initial parameters:")
        print(f"  sp_power={x0[0]:.3f} kW/m³, sp_gas_rate={x0[1]:.3f} cm/s, frother_conc={x0[2]:.1f} mg/L")
        print(f"  ret_time={x0[3]:.1f} min, num_cells={int(x0[4])}")
        print(f"  air_fraction={x0[5]:.3f}, slurry_fraction={x0[6]:.3f}")
        print(f"  bubble_z_pot={x0[7]:.3f} V, cell_volume={x0[8]:.1f} m³, froth_height={x0[9]:.3f} m")

        initial_recovery = -self.objective_function(x0)
        print(f"Initial copper recovery: {initial_recovery*100:.2f}%\n")

        if method == 'differential_evolution':
            # Global optimization using differential evolution
            result = differential_evolution(
                self.objective_function,
                bounds,
                maxiter=100,
                popsize=15,
                tol=0.001,
                atol=0.0001,
                seed=42,
                workers=1,
                disp=True
            )
        elif method == 'slsqp':
            # Local optimization using SLSQP
            result = minimize(
                self.objective_function,
                x0,
                method='SLSQP',
                bounds=bounds,
                options={'maxiter': 100, 'disp': True}
            )
        elif method == 'trust-constr':
            # Local optimization using trust-constr
            result = minimize(
                self.objective_function,
                x0,
                method='trust-constr',
                bounds=bounds,
                options={'maxiter': 100, 'disp': True}
            )
        else:
            raise ValueError(f"Unknown optimization method: {method}")

        # Extract optimal parameters
        optimal_params = {
            'sp_power': result.x[0],
            'sp_gas_rate': result.x[1],
            'frother_conc': result.x[2],
            'ret_time': result.x[3],
            'num_cells': int(round(result.x[4])),
            'air_fraction': result.x[5],
            'slurry_fraction': result.x[6],
            'bubble_z_pot': result.x[7],
            'cell_volume': result.x[8],
            'froth_height': result.x[9],
            'optimal_recovery': -result.fun,
            'success': result.success,
            'message': result.message if hasattr(result, 'message') else 'Optimization complete'
        }

        return optimal_params

    def sensitivity_analysis(self, optimal_params, param_name, param_range):
        """
        Perform sensitivity analysis on a single parameter.

        Args:
            optimal_params: Dict of optimal parameters
            param_name: Name of parameter to vary
            param_range: Array of values to test

        Returns:
            tuple: (param_values, recovery_values)
        """
        param_map = {
            'sp_power': 0,
            'sp_gas_rate': 1,
            'frother_conc': 2,
            'ret_time': 3,
            'num_cells': 4,
            'air_fraction': 5,
            'slurry_fraction': 6,
            'bubble_z_pot': 7,
            'cell_volume': 8,
            'froth_height': 9
        }

        param_idx = param_map[param_name]
        x_optimal = [
            optimal_params['sp_power'],
            optimal_params['sp_gas_rate'],
            optimal_params['frother_conc'],
            optimal_params['ret_time'],
            optimal_params['num_cells'],
            optimal_params['air_fraction'],
            optimal_params['slurry_fraction'],
            optimal_params['bubble_z_pot'],
            optimal_params['cell_volume'],
            optimal_params['froth_height']
        ]

        recovery_values = []
        for value in param_range:
            x_test = x_optimal.copy()
            x_test[param_idx] = value
            recovery = -self.objective_function(x_test)
            recovery_values.append(recovery)

        return param_range, recovery_values


def main():
    """
    Main function to demonstrate optimization.
    """
    print("="*100)
    print("FLOTATION RECOVERY OPTIMIZATION MODEL")
    print("="*100)

    # Define feed composition (from physical_parameters.py)
    feed_composition = {
        'size_classes': [212.13, 106.07, 38.73, 10.00],
        'grade_classes': [0.05, 0.20, 0.40, 0.75, 1.0],
        'throughput_distribution': [
            [136.22, .87, .42, .20, .29],
            [156.95, .61, .18, .14, .12],
            [309.24, .92, .61, .38, 1.85],
            [386.63, .76, .34, .39, 3.88],
        ],
        'sg_by_grade': [3.6, 4, 4.2, 4.2, 4.5],
        'contact_angles_by_grade': [10, 12.5, 27, 52.3, 60],
        'particle_z_pot_by_grade': [-50, -50, -50, -50, -50],
        'grade_by_grade': [0.002, 8, 12, 28, 34.65]
    }

    # Constant parameters (not optimized)
    constant_params = {
        'frother': 4,  # Octanol (2=MIBC, 3=PPG400, 4=Octanol, 5=Pentanol)
        'permitivity': 8.854e-12,  # F/m
        'dielectric': 80,  # Dielectric constant (78-80 for water)
        'pe': 4,  # Peclet number
        'water_or_particle': "Particle",
        'dbl_cell_area': 200,  # Cell area (m^2)
        'dbl_bbl_ratio': 10  # Bubble ratio
    }

    # Create optimizer
    optimizer = FlotationOptimizer(feed_composition, constant_params)

    # Run optimization
    optimal_params = optimizer.optimize(method='differential_evolution')

    # Display results
    print("\n" + "="*100)
    print("OPTIMIZATION RESULTS")
    print("="*100)
    print(f"Optimization successful: {optimal_params['success']}")
    print(f"Message: {optimal_params['message']}")
    print(f"\nOptimal Parameters:")
    print(f"  Specific Power:        {optimal_params['sp_power']:.3f} kW/m³")
    print(f"  Specific Gas Rate:     {optimal_params['sp_gas_rate']:.3f} cm/s")
    print(f"  Frother Concentration: {optimal_params['frother_conc']:.2f} mg/L")
    print(f"  Retention Time:        {optimal_params['ret_time']:.2f} minutes")
    print(f"  Number of Cells:       {optimal_params['num_cells']}")
    print(f"  Air Fraction:          {optimal_params['air_fraction']:.3f} ({optimal_params['air_fraction']*100:.1f}%)")
    print(f"  Slurry Fraction:       {optimal_params['slurry_fraction']:.3f} ({optimal_params['slurry_fraction']*100:.1f}%)")
    print(f"  Bubble Zeta Potential: {optimal_params['bubble_z_pot']:.4f} V")
    print(f"  Cell Volume:           {optimal_params['cell_volume']:.1f} m³")
    print(f"  Froth Height:          {optimal_params['froth_height']:.3f} m ({optimal_params['froth_height']*100:.1f} cm)")
    print(f"\nOptimal Copper Recovery: {optimal_params['optimal_recovery']*100:.2f}%")
    print("="*100)

    # Sensitivity analysis example
    print("\nPerforming sensitivity analysis on specific power...")
    sp_power_range = np.linspace(0.8, 1.3, 20)
    param_values, recovery_values = optimizer.sensitivity_analysis(
        optimal_params,
        'sp_power',
        sp_power_range
    )

    print(f"\n{'Sp. Power (kW/m³)':<20} {'Cu Recovery (%)':<15}")
    print("-"*35)
    for p, r in zip(param_values[::4], recovery_values[::4]):  # Show every 4th value
        print(f"{p:<20.3f} {r*100:<15.2f}")

    return optimal_params


if __name__ == "__main__":
    optimal_params = main()