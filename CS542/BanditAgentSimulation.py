#! usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

class Bandit:
    def __init__(self, k, means, sdev):
        self.k = k
        self.arms = self.generate_arms(k, means, sdev)

    def generate_arms(self, k, mean_list, sdev_list):
        # Each arm contains a tuple containing its mean and standard deviation.
        arms = {}
        for i in range(k):
            arms[i+1] = (mean_list[i], sdev_list[i])
        return arms
    
    def pull(self, i):
        return np.random.normal(loc=self.arms[i][0], scale=self.arms[i][1])
    
class Agent:
    def __init__(self, arms):
        # Track both expected value and number of times arm was selected.
        self.expected_value = self.initialize_expected_value(arms)
        self.action_track = self.initialize_action_track(arms)

    def initialize_expected_value(self, bandit_arms):
        return {key:0 for key in bandit_arms.keys()}
    
    def initialize_action_track(self, bandit_arms):
        return {key:0 for key in bandit_arms.keys()}

    def update_action_value(self, i, r):
        # Increment the action count
        self.action_track[i] += 1
        
        # Update the expected value using incremental averaging
        self.expected_value[i] += (r - self.expected_value[i]) / self.action_track[i]

    def choose_action(self):
        # Choose from the arms with the highest expected value.
        best_arms = [key for key, value in self.expected_value.items() if value == max(self.expected_value.values())]
        return np.random.choice(best_arms)
    
def main():
    # Initialize bandit and agent.
    k = 10
    means = np.random.uniform(-0.5, 0.5, k)
    sdevs = np.ones(k)

    bandit = Bandit(k, means, sdevs)
    bandit.arms

    agent = Agent(bandit.arms)
    agent.expected_value

    # Initialize number of time steps and rewards.
    timesteps = 1000
    reward_at_timestep = np.zeros(timesteps)

    for t in range(timesteps):
        # Select arm to pull.
        i = agent.choose_action()

        # Update action value according to reward received.
        r = bandit.pull(i)
        agent.update_action_value(i, r)

        # Record reward.
        reward_at_timestep[t] = r

    # Calculate average reward per time step.
    cumulative_rewards = np.cumsum(reward_at_timestep)
    average_reward_per_timestep = cumulative_rewards / np.arange(1, timesteps + 1)

    plt.plot(np.arange(1, timesteps + 1), average_reward_per_timestep)
    plt.xlabel('Time Step')
    plt.ylabel('Average Reward')
    plt.title('Average Reward Per Time Step Using Greedy Strategy')
    plt.savefig('BanditPlot.png')
    plt.show()

if __name__ == "__main__":
    main()