import numpy as np

# region - helper functions
def sigmoid(x, b=1, a=0):
    x = np.asarray(x)
    return 1 / (1 + np.exp(-(a + b * x)))

# endregion

# region - models
def Case_1(subj_data, parameters):
    """
    Simulate behavior in a two-armed bandit task with full feedback,
    for a diminishing learning-rate model with increasing choice determinism.

    Expected values are updated using a learning rate that decreases as
    observations accumulate within a block. Choice probabilities are
    generated using a sigmoid choice rule whose inverse temperature
    increases linearly with the number of observations.

    Parameters
    ----------
    subj_data : list[pandas.DataFrame]
        Input dataset. The first element must contain:

        - trial : trial number (used to reset learning at the beginning of each block)
        - b1    : outcome associated with option 1
        - b2    : outcome associated with option 2

    parameters : dict
        Model parameters:

        - beta_change : rate of increase in inverse temperature

    Returns
    -------
    simulated_data : list[pandas.DataFrame]
        Copy of the input data containing:

        - sim_choice_prob : simulated probability of choosing option 1
        - Q1             : value estimate for option 1
        - Q2             : value estimate for option 2

    Notes
    -----
    Used as the generative model in Case 1.
    """

    # Unpack data
    trial = subj_data[0]['trial']
    b1 = subj_data[0]['b1']
    b2 = subj_data[0]['b2']

    n_rows = len(subj_data[0])

    # Model parameters
    prior = 0.5
    beta_change = parameters['beta_change']

    # Initialize arrays
    Q1 = np.zeros(n_rows + 1) * np.nan
    Q2 = np.zeros(n_rows + 1) * np.nan
    sim_choices = np.zeros(n_rows)

    N = 1

    # Simulate learning and choice
    for i in range(n_rows):

        # Reset values at the start of each block
        if trial[i] == 1:
            Q1[i], Q2[i] = prior, prior
            N = 1

        # Choice sensitivity increases with experience
        sim_choices[i] = sigmoid(
            Q1[i] - Q2[i],
            beta_change * N
        )

        # Running-average update rule
        lr_i = 1 / N

        Q1[i + 1] = Q1[i] + lr_i * (b1[i] - Q1[i])
        Q2[i + 1] = Q2[i] + lr_i * (b2[i] - Q2[i])

        N += 1

    # Store simulation output
    simulated_data = subj_data.copy()

    simulated_data[0]['sim_choice_prob'] = sim_choices
    simulated_data[0]['Q1'] = Q1[:n_rows]
    simulated_data[0]['Q2'] = Q2[:n_rows]

    return simulated_data

def Case_2(subj_data, parameters):
    """
    Simulate behavior in a two-armed bandit task with partial feedback,
    for a diminishing learning rate model with fixed choice determinism.

    Expected values are updated only for the chosen option using a
    learning rate that decreases as observations of that option
    accumulate. Choice probabilities are generated using a sigmoid
    choice rule with a fixed inverse temperature.

    Parameters
    ----------
    subj_data : list[pandas.DataFrame]
        Input dataset. The first element must contain:

        - trial : trial number (used to reset learning at the beginning of each block)
        - b1    : outcome associated with option 1
        - b2    : outcome associated with option 2

    parameters : dict
        Model parameters:

        - beta : fixed inverse temperature

    Returns
    -------
    simulated_data : list[pandas.DataFrame]
        Copy of the input data containing:

        - sim_choice_prob : simulated probability of choosing option 1
        - sim_choice      : simulated binary choice
        - Q1              : value estimate for option 1
        - Q2              : value estimate for option 2

    Notes
    -----
    Used as the generative model in Case 2. Learning rates decrease
    separately for each option according to the number of times that
    option has been chosen.
    """

    # Unpack data
    trial = subj_data[0]['trial']
    b1 = subj_data[0]['b1']
    b2 = subj_data[0]['b2']

    n_rows = len(subj_data[0])

    # Model parameters
    prior = 0.5
    beta = parameters['beta']

    # Initialize arrays
    Q1 = np.zeros(n_rows + 1) * np.nan
    Q2 = np.zeros(n_rows + 1) * np.nan

    sim_choices = np.zeros(n_rows) * np.nan
    choices = np.zeros(n_rows) * np.nan

    # Separate observation counts for each option
    N = [1, 1]

    # Simulate learning and choice
    for i in range(n_rows):

        # Reset values at the start of each block
        if trial[i] == 1:
            Q1[i], Q2[i] = prior, prior
            N = [1, 1]

        # Choice probability under fixed inverse temperature
        sim_choices[i] = sigmoid(Q1[i] - Q2[i], beta)

        # Sample a choice
        choices[i] = sim_choices[i] > np.random.uniform(0, 1)

        if choices[i] == 1:

            # Update chosen option only
            lr = 1 / N[0]

            Q1[i + 1] = Q1[i] + lr * (b1[i] - Q1[i])
            Q2[i + 1] = Q2[i]

            N[0] += 1

        else:

            # Update chosen option only
            lr = 1 / N[1]

            Q1[i + 1] = Q1[i]
            Q2[i + 1] = Q2[i] + lr * (b2[i] - Q2[i])

            N[1] += 1

    # Store simulation output
    simulated_data = subj_data.copy()

    simulated_data[0]['sim_choice_prob'] = sim_choices
    simulated_data[0]['sim_choice'] = choices
    simulated_data[0]['Q1'] = Q1[:n_rows]
    simulated_data[0]['Q2'] = Q2[:n_rows]

    return simulated_data

def Case_2_extra(subj_data, parameters):
    """
    Simulate behavior in a two-armed bandit task with full feedback,
    diminishing learning rates, and fixed choice determinism.

    Expected values for both options are updated on every trial using a
    learning rate that decreases as the number of observations increases.
    In contrast to the partial-feedback version (Case 2), both option values
    are updated on every trial since outcomes for both options are observed.

    Choice probabilities are generated using a sigmoid choice rule applied
    to the difference between expected values, with a fixed inverse temperature.

    This model serves as a full-information counterpart to Case 2, designed
    to demonstrate that the directional effects observed in Case 2 arise from
    the partial-feedback structure rather than from the learning rule itself.

    Parameters
    ----------
    subj_data : list[pandas.DataFrame]
        Input dataset. The first element must contain:

        - trial : trial number (used to reset learning at the beginning of each block)
        - b1    : outcome associated with option 1
        - b2    : outcome associated with option 2

    parameters : dict
        Model parameters:

        - beta : fixed inverse temperature

    Returns
    -------
    simulated_data : list[pandas.DataFrame]
        Copy of the input data containing:

        - sim_choice_prob : simulated probability of choosing option 1
        - Q1             : value estimate for option 1
        - Q2             : value estimate for option 2

    Note
    -----
    This is a full-information extension of Case 2. It uses the same learning
    rule but updates both option values on every trial. This extension
    demonstrates that the direction of the systematically positive bias observed in Case 2
    arises from the partial-information structure.
    """

    # Unpack data
    trial = subj_data[0]['trial']
    b1 = subj_data[0]['b1']
    b2 = subj_data[0]['b2']

    n_rows = len(subj_data[0])

    # Parameters
    prior = 0.5
    beta = parameters['beta']

    # Initialize
    Q1 = np.zeros(n_rows + 1) * np.nan
    Q2 = np.zeros(n_rows + 1) * np.nan

    sim_choices = np.zeros(n_rows)

    n_observations = 1

    # Simulate learning
    for i in range(n_rows):

        # Reset at block start
        if trial[i] == 1:
            Q1[i], Q2[i] = prior, prior
            n_observations = 1

        # Choice probability (fixed inverse temperature)
        sim_choices[i] = sigmoid(Q1[i] - Q2[i], beta)

        # Diminishing learning rate
        lr_i = 1 / n_observations

        # Full information update: both options updated every trial
        Q1[i + 1] = Q1[i] + lr_i * (b1[i] - Q1[i])
        Q2[i + 1] = Q2[i] + lr_i * (b2[i] - Q2[i])

        n_observations += 1

    # Store output
    simulated_data = subj_data.copy()
    simulated_data[0]['sim_choice_prob'] = sim_choices
    simulated_data[0]['Q1'] = Q1[:n_rows]
    simulated_data[0]['Q2'] = Q2[:n_rows]

    return simulated_data

# endregion

# Model registry for easy access in data simulation
model_data = {'Case_1':          {'model': Case_1,               'params': ['beta_change']},
              'Case_2':          {'model': Case_2,               'params': ['beta']},
              'Case_2_extra':    {'model': Case_2_extra,         'params': ['beta']}
              }
