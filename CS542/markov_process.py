import numpy as np

def markov_process(num_states, num_steps):
    # Initialize history, states, and starting state.
    history = []
    states = np.arange(num_states)
    current_state = np.random.choice(states)
    history.append(current_state)

    for _ in range(num_steps):
        # 50/50 chance to say or transition
        choice = np.random.choice(['stay', 'transition'], p=[0.5, 0.5])
        
        if choice == 'stay':
            # Stay in current state.
            history.append(current_state)

        if choice == 'transition':
            # Transition to new state.
            other_states = np.setdiff1d(states, current_state)
            new_state = np.random.choice(other_states)
            current_state = new_state
            history.append(current_state)

    # History list contains initial state and subsequent states for num_steps.
    return history

if __name__ == "__main__":
    print(markov_process(5,10))