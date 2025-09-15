# main.py
from config import get_user_input
from simulation import simulate_heat_equation, simulate_fitzhugh_nagumo, simulate_burgers
from utils import plot_results

def main():
    config = get_user_input()
    print("Configuration in main: ", config.initial_conditions)
    if config.equation_type == 'heat':
        initial_globs = [{'position': pos, 'value': val} for pos, val in config.initial_conditions]
        results = simulate_heat_equation(initial_globs, config)
        print("Simulation Complete! Plotting Results...")
        plot_results(results, config.equation_type)
    elif config.equation_type == 'fitzhugh-nagumo':
        initial_globs = [{'position': pos, 'value': [u, v]} for pos, (u, v) in config.initial_conditions]
        results = simulate_fitzhugh_nagumo(initial_globs, config)
        plot_results(results, config.equation_type)
    elif config.equation_type == 'burgers':
        initial_globs = [{'position': pos, 'value': val} for pos, val in config.initial_conditions]
        print("Initial Globs before simulation: ", initial_globs)
        results = simulate_burgers(initial_globs, config)
        plot_results(results, config.equation_type)



if __name__ == "__main__":
    main()
