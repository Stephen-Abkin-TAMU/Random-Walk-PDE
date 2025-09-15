import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as mcolors

def plot_results(results, equation_type):
    """
    Plot the simulation results for the specified equation type.
    
    :param results: List of dictionaries with position and value keys
    :param equation_type: Type of equation being simulated
    
    :return: None"""
    positions = np.array([glob['position'] for glob in results])
    if equation_type == 'heat':
        values = np.array([glob['value'] for glob in results])
        plt.scatter(positions, values, c=values, cmap='coolwarm', edgecolor='none')
        plt.colorbar(label='Temperature')
        plt.title('Heat Distribution Over Space')
        plt.xlabel('Position')
        plt.ylabel('Temperature')
        plt.grid(True)
        plt.show()

    elif equation_type == 'fitzhugh-nagumo':
        u_values = np.array([glob['value'][0] for glob in results])
        v_values = np.array([glob['value'][1] for glob in results])
        
        plt.figure(figsize=(12, 6))
        ax1 = plt.subplot(1, 2, 1)
        scatter = ax1.scatter(positions, u_values, c=u_values, cmap='coolwarm', edgecolor='none')
        plt.colorbar(scatter, ax=ax1, label='Membrane Potential (u)')
        ax1.set_title('Membrane Potential (u) Over Space')
        ax1.set_xlabel('Position')
        ax1.set_ylabel('Membrane Potential (u)')

        ax2 = plt.subplot(1, 2, 2)
        scatter = ax2.scatter(positions, v_values, c=v_values, cmap='viridis', edgecolor='none')
        plt.colorbar(scatter, ax=ax2, label='Recovery Variable (v)')
        ax2.set_title('Recovery Variable (v) Over Space')
        ax2.set_xlabel('Position')
        ax2.set_ylabel('Recovery Variable (v)')

        plt.suptitle(f'FitzHugh-Nagumo Model Results')
        plt.tight_layout(pad=3.0)

    elif equation_type == 'burgers':
        values = np.array([glob['value'][0] for glob in results])
        plt.figure(figsize=(10, 5))
        plt.plot(positions, values, label='Velocity Profile', color='blue', linewidth=2)
        plt.scatter(positions, values, color='red')  # Points to highlight the data
        plt.title('Burger\'s Equation: Velocity Profile Over Space')
        plt.xlabel('Position')
        plt.ylabel('Velocity')
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
    
    plt.savefig(f'output/{equation_type}_simulation_results.png')
    plt.close()
