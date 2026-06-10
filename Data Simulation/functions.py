import numpy as np
import pandas as pd
from pathlib import Path
from models import model_data


def create_data(model_name, parameters, subject_conditions,
                n_subjects, n_blocks, n_trials,
                save_dir= None, alt_name=None):
    """
    Generate synthetic two-armed bandit data and simulate choices
    using a specified computational model.

    The function proceeds in two stages:
    (1) generate trial-level reward structures for each subject
    (2) simulate choices using a model-defined learning algorithm

    Parameters
    ----------
    model_name : str
        Key identifying the generative model in `model_data`.

    parameters : dict
        Dictionary of model parameters, each entry is a list or array
        of length `n_subjects`.

    subject_conditions : str or list of str
        Condition assignment per subject. If a single string is provided,
        it is replicated across all n_subjects.

    save_dir : str
        Directory where the simulated dataset will be saved.

    n_subjects : int
        Number of n_subjects to simulate.

    n_blocks : int
        Number of n_blocks per block.

    n_blocks : int
        Number of n_blocks per trial block.

    alt_name : str, optional
        Alternative filename for saving output.

    Returns
    -------
    None
        Saves a CSV file containing the full simulated dataset.
    """

    # ------------------------------------------------------------------
    # Ensure subject-level condition assignment
    # ------------------------------------------------------------------
    if isinstance(subject_conditions, str):
        subject_conditions = [subject_conditions] * n_subjects

    # ------------------------------------------------------------------
    # Generate reward structure (two-armed bandit environment)
    # ------------------------------------------------------------------
    data = pd.DataFrame(
        columns=['subject', 'block', 'trial', 'b1', 'b2', 'choice']
    )

    for s in range(n_subjects):
        for block in range(n_blocks):

            # Define reward probabilities per condition
            if subject_conditions[s] == 'same':
                p = np.random.uniform(0.1, 0.9)
                probs = [p, p]

            elif subject_conditions[s] == 'opposing':
                p = np.random.uniform(0.1, 0.9)
                probs = [p, 1 - p]

            elif subject_conditions[s] == 'constant':
                probs = [0.25, 0.75]

            else:  # random condition
                probs = [
                    np.random.uniform(0.1, 0.9),
                    np.random.uniform(0.1, 0.9)
                ]

            block_data = pd.DataFrame({
                'subject': np.repeat(s + 1, n_trials),
                'condition': np.repeat(subject_conditions[s], n_trials),
                'block': np.repeat(block + 1, n_trials),
                'trial': np.arange(1, n_trials + 1),

                # Bernoulli outcomes for each option
                'b1': np.random.binomial(1, probs[0], n_trials),
                'b2': np.random.binomial(1, probs[1], n_trials),

                # placeholders for model output
                'choice': np.nan,
                'Q1': np.nan,
                'Q2': np.nan
            })

            data = pd.concat([data, block_data], ignore_index=True)

    # Convert to subject-wise dataframes, organized as a list of dataframes (one per subject)
    data_list = [
        [data[data['subject'] == s].reset_index(drop=True)]
        for s in data['subject'].unique()
    ]

    # ------------------------------------------------------------------
    # Simulate choices using selected computational model
    # ------------------------------------------------------------------
    model = model_data[model_name]['model']

    sim_data = pd.DataFrame()

    for s in range(n_subjects):

        subj_sim_data = model(subj_data=data_list[s],
                              parameters={param: parameters[param][s]for param in parameters})

        # Attach subject-specific parameters to df (for traceability)
        for param in parameters:
            if param in model_data[model_name]['params']:
                subj_sim_data[0][param] = parameters[param][s]

        sim_data = pd.concat([sim_data, subj_sim_data[0]], ignore_index=True)

    # ------------------------------------------------------------------
    # Generate final choice variable
    # ------------------------------------------------------------------

    # If the model already outputs a binary choice variable, use it directly
    if 'sim_choice' in sim_data.columns:
        sim_data['choice'] = sim_data['sim_choice']

    # Otherwise, sample choices based on the simulated choice probabilities
    else:
        sim_data['choice'] = np.random.binomial(1,sim_data['sim_choice_prob'],len(sim_data))

    # ------------------------------------------------------------------
    # Final formatting
    # ------------------------------------------------------------------
    sim_data.rename(columns={'sim_choice_prob': 'p_choice',
                             'Q1':              'true_Q1',
                             'Q2':              'true_Q2'},
                    inplace=True)

    # ------------------------------------------------------------------
    # Save output
    # ------------------------------------------------------------------
    if alt_name is not None:
        model_name = alt_name

    if save_dir is None:
        save_dir = Path(__file__).resolve().parent / "Simulated Data"
    else:
        save_dir = Path(save_dir)

    save_dir.mkdir(parents=True, exist_ok=True)

    filename = f"{save_dir}/{model_name}.csv"
    sim_data.to_csv(filename, index=False)

    print(f"Simulated data saved to {filename}")