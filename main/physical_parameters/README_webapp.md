# Copper Production Optimizer - Web Application

A Flask-based web application for optimizing copper extraction processes with physics-based modeling.

## Features

- **Interactive Input Sliders**
  - Feed ore grade: 0.2% to 1%
  - Target copper production: 1,000 to 100,000 tons/year

- **Optimization Results**
  - NPV comparison (with and without electrolyzer)
  - IRR and payback period
  - Optimal operating parameters for all process stages
  - Recovery rates for each stage

## Local Development

### Install Dependencies

```bash
pip install -r requirements.txt
```

### Run the Application

```bash
python app.py
```

The application will be available at `http://localhost:5000`

## Deployment to GitHub Pages (Static Hosting)

Note: GitHub Pages only supports static sites. For a Flask app, you'll need to use a service that supports Python applications.

### Option 1: Deploy to Heroku (Recommended)

1. Create a `Procfile`:
```
web: gunicorn app:app
```

2. Add `gunicorn` to requirements.txt

3. Deploy to Heroku:
```bash
heroku create copper-optimizer
git push heroku main
```

### Option 2: Deploy to Render

1. Create account at render.com
2. Connect your GitHub repository
3. Set build command: `pip install -r requirements.txt`
4. Set start command: `gunicorn app:app`

### Option 3: Deploy to PythonAnywhere

1. Create account at pythonanywhere.com
2. Upload your files
3. Configure WSGI file to point to your Flask app
4. Set working directory

## File Structure

```
physical_parameters/
├── app.py                 # Flask application
├── optimal_tea.py         # Optimization engine
├── process_cost_integration.py
├── flotation_recovery.py
├── leaching.py
├── solvent_extraction.py
├── tea.py
├── templates/
│   └── index.html        # Web interface
├── requirements.txt      # Python dependencies
└── README_webapp.md      # This file
```

## API Endpoint

### POST /optimize

Request body:
```json
{
  "feed_grade": 0.0058,
  "target_tonnage": 10000,
  "max_iterations": 30
}
```

Response:
```json
{
  "success": true,
  "without_electrolyzer": {
    "npv": 123456789,
    "irr": 15.5,
    "payback_years": 6.2,
    ...
  },
  "with_electrolyzer": {
    "npv": 145678901,
    "irr": 18.2,
    "payback_years": 5.1,
    ...
  },
  "comparison": {
    "npv_improvement": 22222112,
    "npv_improvement_pct": 18.0,
    "recommendation": "Use Electrolyzer"
  }
}
```

## Notes

- Optimization typically takes 30-60 seconds depending on parameters
- The app uses differential evolution optimization with 30 iterations by default
- All costs are in USD, copper price is fixed at $13,000/ton
- Concentrate grade is fixed at 28% (industry standard)
