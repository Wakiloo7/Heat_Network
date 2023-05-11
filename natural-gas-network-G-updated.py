import pyomo.environ as pyo

# Define the sets
phi_k = {1, 2, 3} #Represents time periods
phi_n = {1, 2, 3} #Represents gas nodes
phi_gmn = {(1, 2), (2, 3)} #Represents gas edges, which are pairs of connected nodes

# Define the parameters
C_n = {1: 1.0, 2: 1.2, 3: 1.8} # Represents the cost coefficient at node n
P_n_min = {1: 10, 2: 15, 3: 20} # Minimum gas pressure at node n
P_n_max = {1: 100, 2: 90, 3: 80} # Maximum gas pressure at node n

G_n_min = {1: 0, 2: 0, 3: -100}  # Minimum gas flow at node n
G_n_max = {1: 100, 2: 100, 3: 0} # Maximum gas flow at node n

gf_mn_max = {(1, 2): 1000, (2, 3): 800} #Maximum gas flow between nodes m and n
C_mn=18


# Create the model
model = pyo.ConcreteModel()

# Define the variables
model.gf_kmn = pyo.Var(phi_k, phi_gmn, bounds=lambda model, k, m, n: (-gf_mn_max[m, n], gf_mn_max[m, n]))
model.P_kn = pyo.Var(phi_k, phi_n, initialize=lambda model, k, n: P_n_min[n])
model.G_kn = pyo.Var(phi_k, phi_n, initialize=lambda model, k, n: max(G_n_min[n], 0))

# Define the binary variables I_kmn_plus and I_kmn_minus
model.I_kmn_plus = pyo.Var(phi_k, phi_gmn, within=pyo.Binary)
model.I_kmn_minus = pyo.Var(phi_k, phi_gmn, within=pyo.Binary)

# Add constraints for the relationship between I_kmn_plus, I_kmn_minus, and gf_kmn
def gas_flow_direction_positive(model, k, m, n):
    return model.gf_kmn[k, m, n] <= model.I_kmn_plus[k, m, n] * gf_mn_max[m, n]

def gas_flow_direction_negative(model, k, m, n):
    return model.gf_kmn[k, m, n] >= model.I_kmn_minus[k, m, n] * (-gf_mn_max[m, n])


def gas_flow_direction_sum(model, k, m, n):
    return model.I_kmn_plus[k, m, n] + model.I_kmn_minus[k, m, n] == 1

model.gas_flow_direction_positive_constraint = pyo.Constraint(phi_k, phi_gmn, rule=gas_flow_direction_positive)
model.gas_flow_direction_negative_constraint = pyo.Constraint(phi_k, phi_gmn, rule=gas_flow_direction_negative)
model.gas_flow_direction_sum_constraint = pyo.Constraint(phi_k, phi_gmn, rule=gas_flow_direction_sum)


#Define the gas volume balance constraint
def gas_volume_balance_rule(model, k, n):
    inflow = sum(model.gf_kmn[k, m, n] for m in phi_n if (m, n) in phi_gmn)
    outflow = sum(model.gf_kmn[k, n, m] for m in phi_n if (n, m) in phi_gmn)
    return inflow - outflow - model.G_kn[k, n] == 0

model.gas_volume_balance_constraint = pyo.Constraint(phi_k, phi_n, rule=gas_volume_balance_rule)


#Define auxiliary variable Gamma_kmn
model.Gamma_kmn = pyo.Var(phi_k, phi_gmn, within=pyo.NonNegativeReals)

#McCormick envelope constraints
def mccormick1(model, k, m, n):
    return model.Gamma_kmn[k, m, n] >= (1 / C_mn) ** 2 * model.gf_kmn[k, m, n] ** 2

def mccormick2(model, k, m, n):
    return model.Gamma_kmn[k, m, n] >= model.P_kn[k, n] - model.P_kn[k, m] + (model.I_kmn_plus[k, m, n] - model.I_kmn_minus[k, m, n] + 1) * (P_n_min[m] - P_n_max[n])

def mccormick3(model, k, m, n):
    return model.Gamma_kmn[k, m, n] >= model.P_kn[k, m] - model.P_kn[k, n] + (model.I_kmn_plus[k, m, n] - model.I_kmn_minus[k, m, n] - 1) * (P_n_max[m] - P_n_min[n])

def mccormick4(model, k, m, n):
    return model.Gamma_kmn[k, m, n] <= model.P_kn[k, n] - model.P_kn[k, m] + (model.I_kmn_plus[k, m, n] - model.I_kmn_minus[k, m, n] + 1) * (P_n_max[m] - P_n_min[n])

def mccormick5(model, k, m, n):
    return model.Gamma_kmn[k, m, n] <= model.P_kn[k, m] - model.P_kn[k, n] + (model.I_kmn_plus[k, m, n] - model.I_kmn_minus[k, m, n] - 1) * (P_n_min[m] - P_n_max[n])

model.mccormick1_constraint = pyo.Constraint(phi_k, phi_gmn, rule=mccormick1)
model.mccormick2_constraint = pyo.Constraint(phi_k, phi_gmn, rule=mccormick2)
model.mccormick3_constraint = pyo.Constraint(phi_k, phi_gmn, rule=mccormick3)
model.mccormick4_constraint = pyo.Constraint(phi_k, phi_gmn, rule=mccormick4)
model.mccormick5_constraint = pyo.Constraint(phi_k, phi_gmn, rule=mccormick5)

#Add constraints for minimum and maximum pressure values
def pressure_min_constraint(model, k, n):
    return model.P_kn[k, n] >= P_n_min[n]

def pressure_max_constraint(model, k, n):
    return model.P_kn[k, n] <= P_n_max[n]

model.pressure_min_constraint = pyo.Constraint(phi_k, phi_n, rule=pressure_min_constraint)
model.pressure_max_constraint = pyo.Constraint(phi_k, phi_n, rule=pressure_max_constraint)

# Define constraints for minimum and maximum gas values
# def gas_min_constraint(model, k, n):
#     if n != 3:
#         return model.G_kn[k, n] >= G_n_min[n]
#     return pyo.Constraint.Skip

# def gas_max_constraint(model, k, n):
#     if n != 3:
#         return model.G_kn[k, n] <= G_n_max[n]
#     return pyo.Constraint.Skip

# # Add constraints to the model
# model.gas_min_constraint = pyo.Constraint(phi_k, phi_n, rule=gas_min_constraint)
# model.gas_max_constraint = pyo.Constraint(phi_k, phi_n, rule=gas_max_constraint)

# Define constraints for minimum and maximum gas values
def gas_min_constraint(model, k, n):
    return model.G_kn[k, n] >= G_n_min[n]

def gas_max_constraint(model, k, n):
    return model.G_kn[k, n] <= G_n_max[n]

model.gas_min_constraint = pyo.Constraint(phi_k, phi_n, rule=gas_min_constraint)
model.gas_max_constraint = pyo.Constraint(phi_k, phi_n, rule=gas_max_constraint)


# Define the objective function
def objective_function(model):
    return sum(C_n[n] * model.G_kn[k, n] for k in phi_k for n in phi_n)

model.objective = pyo.Objective(rule=objective_function, sense=pyo.minimize)

#Solve the model
solver = pyo.SolverFactory('gurobi')
solver.options['MIPGap'] = 1e-6  # Set the MIP gap to a smaller value (default is 1e-4)
results = solver.solve(model, tee=True, options={'iis': 'find'})


# Check if the solver terminated with an optimal solution
if results.solver.status == pyo.SolverStatus.ok and results.solver.termination_condition == pyo.TerminationCondition.optimal:
    print("Solver terminated with an optimal solution")
    # Print the results
    print("\nGas Flows (gf_kmn):")
    for k in phi_k:
        for m, n in phi_gmn:
            print(f"k: {k}, m: {m}, n: {n}, value: {model.gf_kmn[k, m, n].value}")
    print("\nPressures (P_kn):")
    for k in phi_k:
        for n in phi_n:
            print(f"k: {k}, n: {n}, value: {model.P_kn[k, n].value}")

    print("\nGas Volumes (G_kn):")
    for k in phi_k:
        for n in phi_n:
            print(f"k: {k}, n: {n}, value: {model.G_kn[k, n].value}")

    print("\nObjective Function Value:")
    print(model.objective())
else:
    print("Solver did not terminate with an optimal solution.")


