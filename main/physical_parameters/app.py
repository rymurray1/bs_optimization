"""
Flask web application for Copper Production Optimization
Provides a simple interface to optimize process parameters and compare electrolyzer economics
"""

from flask import Flask, render_template, request, jsonify
import sys
import os

# Add current directory to path
sys.path.insert(0, os.path.dirname(__file__))

from optimal_tea import OptimalTEA, calculate_npv, calculate_irr, calculate_payback_period

app = Flask(__name__)

@app.route('/')
def index():
    """Render the main page with input sliders"""
    return render_template('index.html')

@app.route('/optimize', methods=['POST'])
def optimize():
    """
    Run optimization with user inputs and return results.

    Expects JSON with:
        - feed_grade: float (0.002 to 0.01, representing 0.2% to 1%)
        - target_tonnage: int (copper production in tons/year)

    Returns JSON with:
        - Without electrolyzer: NPV, IRR, payback, parameters
        - With electrolyzer: NPV, IRR, payback, parameters
        - NPV improvement
    """
    try:
        data = request.get_json()

        # Extract inputs
        feed_grade = float(data.get('feed_grade', 0.0058))  # Default 0.58%
        target_tonnage = int(data.get('target_tonnage', 10000))  # Default 10,000 tons/year
        max_iterations = int(data.get('max_iterations', 30))  # Optimization iterations

        # Fixed parameters
        concentrate_grade = 0.28  # 28% concentrate grade (industry standard)

        # Run optimization WITHOUT electrolyzer
        optimizer_no = OptimalTEA(
            target_cu_tons=target_tonnage,
            feed_grade=feed_grade,
            concentrate_grade=concentrate_grade,
            use_electrolyzer=False
        )
        results_no = optimizer_no.run_optimization(method='differential_evolution', maxiter=max_iterations)
        result_no = results_no['result']

        # Calculate financial metrics
        npv_no = calculate_npv(result_no['total_capex'], result_no['total_opex'], result_no['cu_produced_tons'])
        irr_no = calculate_irr(result_no['total_capex'], result_no['total_opex'], result_no['cu_produced_tons'])
        payback_no = calculate_payback_period(result_no['total_capex'], result_no['total_opex'], result_no['cu_produced_tons'])

        # Run optimization WITH electrolyzer
        optimizer_yes = OptimalTEA(
            target_cu_tons=target_tonnage,
            feed_grade=feed_grade,
            concentrate_grade=concentrate_grade,
            use_electrolyzer=True
        )
        results_yes = optimizer_yes.run_optimization(method='differential_evolution', maxiter=max_iterations)
        result_yes = results_yes['result']

        # Calculate financial metrics
        npv_yes = calculate_npv(result_yes['total_capex'], result_yes['total_opex'], result_yes['cu_produced_tons'])
        irr_yes = calculate_irr(result_yes['total_capex'], result_yes['total_opex'], result_yes['cu_produced_tons'])
        payback_yes = calculate_payback_period(result_yes['total_capex'], result_yes['total_opex'], result_yes['cu_produced_tons'])

        # Calculate improvement
        npv_improvement = npv_yes - npv_no
        npv_improvement_pct = (npv_improvement / npv_no * 100) if npv_no != 0 else 0

        # Unpack optimal parameters
        grinding_no, flotation_no, leaching_no, sx_no, ew_no = optimizer_no.unpack_parameters(results_no['optimal_params'])
        grinding_yes, flotation_yes, leaching_yes, sx_yes, ew_yes = optimizer_yes.unpack_parameters(results_yes['optimal_params'])

        # Build response
        response = {
            'success': True,
            'inputs': {
                'feed_grade_pct': feed_grade * 100,
                'target_tonnage': target_tonnage,
                'concentrate_grade_pct': concentrate_grade * 100
            },
            'without_electrolyzer': {
                'npv': npv_no,
                'irr': irr_no * 100,  # Convert to percentage
                'payback_years': payback_no,
                'total_capex': result_no['total_capex'],
                'total_opex': result_no['total_opex'],
                'copper_produced': result_no['cu_produced_tons'],
                'ore_feed_required': result_no['ore_feed_tpa'],
                'recoveries': {
                    'flotation': result_no['recoveries']['flotation'] * 100,
                    'leaching': result_no['recoveries']['leaching'] * 100,
                    'sx': result_no['recoveries']['sx'] * 100,
                    'ew': result_no['recoveries']['ew'] * 100,
                    'overall': result_no['recoveries']['overall'] * 100
                },
                'optimal_parameters': {
                    'grinding': grinding_no,
                    'flotation': flotation_no,
                    'leaching': leaching_no,
                    'sx': sx_no,
                    'ew': ew_no
                }
            },
            'with_electrolyzer': {
                'npv': npv_yes,
                'irr': irr_yes * 100,
                'payback_years': payback_yes,
                'total_capex': result_yes['total_capex'],
                'total_opex': result_yes['total_opex'],
                'copper_produced': result_yes['cu_produced_tons'],
                'ore_feed_required': result_yes['ore_feed_tpa'],
                'recoveries': {
                    'flotation': result_yes['recoveries']['flotation'] * 100,
                    'leaching': result_yes['recoveries']['leaching'] * 100,
                    'sx': result_yes['recoveries']['sx'] * 100,
                    'ew': result_yes['recoveries']['ew'] * 100,
                    'overall': result_yes['recoveries']['overall'] * 100
                },
                'optimal_parameters': {
                    'grinding': grinding_yes,
                    'flotation': flotation_yes,
                    'leaching': leaching_yes,
                    'sx': sx_yes,
                    'ew': ew_yes
                }
            },
            'comparison': {
                'npv_improvement': npv_improvement,
                'npv_improvement_pct': npv_improvement_pct,
                'recommendation': 'Use Electrolyzer' if npv_improvement > 0 else 'Buy Acid'
            }
        }

        return jsonify(response)

    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        }), 500

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)
