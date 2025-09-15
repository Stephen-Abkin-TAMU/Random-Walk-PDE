import numpy as np

class SimulationConfig:
    """
    A class to hold the configuration for the simulation.

    Attributes:
        equation_type (str): The type of PDE to simulate (Heat, Fitzhugh-Nagumo, Burgers).
        domain_type (str): The nature of the domain (Finite, Semi-Infinite, Infinite).
        domain_size (float): The size of the domain.
        boundary_conditions (dict): Conditions at the boundaries of the domain.
        diff_constant (float): Diffusivity constant used in the simulations.
        time_step (float): The time increment for the simulation.
        total_time (float): Total time for which the simulation runs.
        num_points (int): The number of spatial points (globs) in the simulation.
        initial_conditions (list): Initial conditions for the simulation.
        reaction_term (bool): Whether a reaction term is included in the simulation.
        a, b, tau (float, optional): Parameters specific to the Fitzhugh-Nagumo model.
    """
    def __init__(self, equation_type, domain_type, domain_size, boundary_conditions, diff_constant, time_step, total_time, num_points, initial_conditions, reaction_term, a=None, b=None, tau=None):
        self.equation_type = equation_type
        self.domain_type = domain_type
        self.domain_size = domain_size
        self.boundary_conditions = boundary_conditions
        self.diff_constant = diff_constant
        self.time_step = time_step
        self.total_time = total_time
        self.num_points = num_points
        self.initial_conditions = initial_conditions
        self.reaction_term = reaction_term
        self.a = a  # FitzHugh-Nagumo parameter
        self.b = b  # FitzHugh-Nagumo parameter
        self.tau = tau  # FitzHugh-Nagumo temporal scaling parameter




def generate_heat_equation_initial_conditions(domain_size, num_points):
    """
    Generate intial conditions for the heat equation.

    :param domain_size: Size of the spatial domain
    :param num_points: Number of spatial points (globs)
    :return: List of tuples with positions and initial values
    """
    positions = np.linspace(0, domain_size, num_points)
    values = np.linspace(0, 100, num_points) # temperature from 0 to 100
    return list(zip(positions, values))

def generate_fitzhugh_nagumo_initial_conditions(domain_size, num_points, stimulated_region=None, stimulus_magnitude=2):
    """
    Generate initial conditions for FitzHugh-Nagumo model, with an optional stimulus.
    

    :param domain_size: Size of the spatial domain
    :param num_points: Number of spatial points (globs)
    :param stimulated_region: Tuple indicating the start and end of the stimulated region (e.g., (45, 55))
    :param stimulus_magnitude: The initial value for u in the stimulated region
    :return: List of tuples with positions and initial values (u, v)
    """
    positions = np.linspace(0, domain_size, num_points)
    # Initially, all cells are at rest with u = 0 and v = 0
    values = [(0, 0) for _ in positions]

    if stimulated_region:
        start_index = int((stimulated_region[0] / domain_size) * num_points)
        end_index = int((stimulated_region[1] / domain_size) * num_points)
        for i in range(start_index, end_index):
            values[i] = (stimulus_magnitude, 0)  # Apply stimulus magnitude to u, v remains 0

    return list(zip(positions, values))


def generate_burgers_initial_conditions(domain_size, num_points, condition_type='shock', shock_position=None, left_value=1, right_value=0):
    """
    Generate initial conditions for the Burgers' equation.

    :param domain_size: Size of the spatial domain
    :param num_points: Number of spatial points (globs)
    :param condition_type: Type of initial condition ('shock', 'step', 'rarefaction', 'uniform', 'tanh')
    :param shock_position: Position in the domain to start the shock (used if 'shock' or 'step')
    :param left_value: Value before the shock or step position
    :param right_value: Value after the shock or step position
    :return: List of tuples with positions and values
    """
    positions = np.linspace(0, domain_size, num_points)
    values = np.zeros(num_points)  # Ensures the array is of numeric type
    
    if condition_type in ['shock', 'step']:
        shock_index = int(shock_position / domain_size * num_points) if shock_position is not None else num_points // 2
        values[:shock_index] = left_value
        values[shock_index:] = right_value
    elif condition_type == 'rarefaction':
        values = np.linspace(left_value, right_value, num_points)
    elif condition_type == 'uniform':
        values.fill(left_value)
    elif condition_type == 'tanh':
        # Example to add a tanh profile for traveling wave scenarios
        x_center = domain_size / 2
        scale = 0.1 * domain_size  # Scale of the tanh slope
        values = 0.5 * (1 - np.tanh((positions - x_center) / scale))
    
    return list(zip(positions, values))  # Returning a list of tuples to maintain consistency


def get_user_input():
    """
    Retrieve user input to configure the simulation.

    :return: SimulationConfig object with the user-defined parameters
    """
    print("Welcome to the Gradient Random Walk Simulation Setup")
    equation_type = input("Choose the equation to simulate (Heat, Fitzhugh-Nagumo, Burgers): ").lower().strip()
    domain_type = input("Enter the domain type (Finite, Semi-Infinite, Infinite): ")
    domain_size = float(input("Enter the domain size: "))

    boundary_conditions = {}
    if domain_type != "Infinite":
        for end in ["LEFT", "RIGHT"]:
            bc_type = input(f"Enter {end} boundary condition type (Dirichlet, Neumann): ")
            bc_value = float(input(f"Enter {end} the boundary value: "))
            boundary_conditions[end] = {'type': bc_type, 'value': bc_value}

    diff_constant = float(input("Enter the diffusivity constant: "))
    time_step = float(input("Enter the time step (delta_t): "))
    total_time = float(input("Enter the total simulation time: "))
    num_points = int(input("Enter the number of globs: "))
    
    initial_conditions, a, b, tau = None, None, None, None
    stimulated_region, stimulus_magnitude = None, None
    if equation_type == "fitzhugh-nagumo":
        a = float(input("Enter FitzHugh-Nagumo parameter a: "))
        b = float(input("Enter FitzHugh-Nagumo parameter b: "))
        tau = float(input("Enter FitzHugh-Nagumo temporal scaling parameter tau: "))
        if input("Apply an initial electrical stimulus? (yes/no): ").lower() == 'yes':
            start = float(input("Enter the start position of the stimulated region (e.g., 45): "))
            end = float(input("Enter the end position of the stimulated region (e.g., 55): "))
            stimulus_magnitude = float(input("Enter the magnitude of the stimulus: "))
            stimulated_region = (start, end)
        initial_conditions = generate_fitzhugh_nagumo_initial_conditions(domain_size, num_points, stimulated_region, stimulus_magnitude)
    elif equation_type == "heat":
        initial_conditions = generate_heat_equation_initial_conditions(domain_size, num_points)
    elif equation_type == "burgers":
        print("Select the type of initial condition:")
        print("1. Shock Wave")
        print("2. Step Function")
        print("3. Rarefaction Wave")
        print("4. Uniform Condition")
        print("5. Hyperbolic Tangent Profile")
        condition_type = input("Enter your choice (1-5): ")
        
        condition_map = {
            '1': 'shock',
            '2': 'step',
            '3': 'rarefaction',
            '4': 'uniform',
            '5': 'tanh'
        }
        
        condition_type = condition_map.get(condition_type, 'shock')  # Default to 'shock' if invalid input
        
        if condition_type in ['shock', 'step']:
            shock_position = float(input("Enter the position for the transition (e.g., 5): "))
            left_value = float(input("Enter the value before the transition: "))
            right_value = float(input("Enter the value after the transition: "))
            initial_conditions = generate_burgers_initial_conditions(domain_size, num_points, condition_type=condition_type, shock_position=shock_position, left_value=left_value, right_value=right_value)
        elif condition_type == 'rarefaction':
            left_value = float(input("Enter the starting value: "))
            right_value = float(input("Enter the ending value: "))
            initial_conditions = generate_burgers_initial_conditions(domain_size, num_points, condition_type='rarefaction', left_value=left_value, right_value=right_value)
        elif condition_type == 'uniform':
            uniform_value = float(input("Enter the uniform value across the domain: "))
            initial_conditions = generate_burgers_initial_conditions(domain_size, num_points, condition_type='uniform', left_value=uniform_value)
        elif condition_type == 'tanh':
            initial_conditions = generate_burgers_initial_conditions(domain_size, num_points, condition_type='tanh')

    else:
        raise ValueError("Invalid equation type")

    reaction_term = input("Is there a reaction term? (yes/no): ").lower() == 'yes'

    return SimulationConfig(equation_type, domain_type, domain_size, boundary_conditions, diff_constant, time_step, total_time, num_points, initial_conditions, reaction_term, a, b, tau)

if __name__ == "__main__":
    config = get_user_input()
    print("Simulation Configuration:")
    print(vars(config))
