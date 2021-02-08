"""
Test BPP with simulated data!
"""
import numpy as np

sample_size = 1765
num_years = 14
years = np.array(range(num_years)) + 1979

simulated_data = np.zeros((num_years*sample_size,9))
for i in range(sample_size):
    this_simulated_data = np.zeros((num_years,9))
    this_simulated_data[:,0] = i
    this_simulated_data[:,1] = years
    this_simulated_data[:,2] = 0
    this_simulated_data[:,3] = 0
    this_simulated_data[8:11,3] = 1
    trans_shocks = 2.0**0.5*np.random.normal(size=num_years+1)
    dif_tran = trans_shocks[1:]-trans_shocks[:-1]
    perm_shocks = np.random.normal(size=num_years)
    this_simulated_data[:,4] = dif_tran + perm_shocks
    phi = 0.8
    psi = 0.3
    this_simulated_data[:,5] = 1
    this_simulated_data[:,6] = phi*perm_shocks + psi*trans_shocks[1:]
    this_simulated_data[:,7] = 1
    this_simulated_data[:,8] = 3
    simulated_data[i*num_years:(i+1)*num_years,:] = this_simulated_data
    
simulated_data = pd.DataFrame(simulated_data)
simulated_data.to_csv(Path("./InputFiles/SimulatedData.csv"), index = False)
    
    
    
