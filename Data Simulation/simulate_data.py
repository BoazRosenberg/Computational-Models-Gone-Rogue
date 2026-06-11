from functions import create_data
from pathlib import Path


if __name__ == "__main__":

    '''
    Main models:
        Case_1, Case_2

    Extra variants of case 2 with full information:
        Case_2_extra_random, Case_2_extra_opposing
    
    '''

    model_names = ['Case_1', 'Case_2', 'Case_2_extra_random', 'Case_2_extra_opposing']

    n_trials = 30
    n_blocks = 30
    n_subjects = 30

    parameters = {'beta':        [3.0] * n_subjects,
                  'beta_change': [0.2] * n_subjects}

    save_dir = Path(__file__).resolve().parent / "Simulated Data"

    for model_name in model_names:

        # Case 2 extra: full-information control (random reward probabilities)
        if model_name == 'Case_2_extra_random':
            create_data('Case_2_extra', parameters, 'random',
                        n_subjects, n_blocks, n_trials, alt_name='Case_2_extra_random')

        # Case 2 extra: full-information control (opposing reward probabilities: [p, 1-p])
        elif model_name == 'Case_2_extra_opposing':
            create_data('Case_2_extra', parameters, 'opposing',
                        n_subjects, n_blocks, n_trials, alt_name='Case_2_extra_opposing')

        else:
            create_data(model_name, parameters, 'random',
                        n_subjects, n_blocks, n_trials)

