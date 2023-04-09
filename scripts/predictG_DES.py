import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import RidgeCV
# Load the cognitive scores data
df = pd.read_csv('/home/users/apines/DES_pcadf_g.csv')

# Path to the rs_BUTD files
rs_BUTD_path = '/oak/stanford/groups/leanew1/users/apines/OpFlAngDs/DES'

# Initialize counter for subjects missing rs_BUTD_L.csv file
missing_files_count = 0

# Create empty arrays to store features and scores
features = []
scores = []

# Loop through each subject and load their features
for idx, row in df.iterrows():
    subject_id = row['id']
    subfolder = f'sub-{subject_id}'
    rs_L_path = f'{rs_BUTD_path}/{subfolder}/'
    rs_R_path = f'{rs_BUTD_path}/{subfolder}/'
    # Check if the rs_BUTD_L and rs_BUTD_R files exist
    try:
        rs_L = pd.read_csv(f'{rs_L_path}rs_BUTD_L.csv')
        rs_L2 = pd.read_csv(f'{rs_L_path}rs2_BUTD_L.csv')
        rs_R = pd.read_csv(f'{rs_R_path}rs_BUTD_R.csv')
        rs_R2 = pd.read_csv(f'{rs_R_path}rs2_BUTD_R.csv')
    except FileNotFoundError:
        # Increment the counter for missing files and continue to the next subject
        missing_files_count += 1
        continue
    # Calculate the average of rs and rs2 separately for L and R
    rs_avg_L = (rs_L['Var1_1'] + rs_L2['Var1_1']) / 2
    rs_avg_R = (rs_R['Var1_1'] + rs_R2['Var1_1']) / 2
    # Combine L and R (averaged over rs and rs2) into one feature vector
    feature = np.concatenate([rs_avg_L, rs_avg_R])
    features.append(feature)
    scores.append(row['nihtb_flanker_fullcorrtscore'])

# Print the number of subjects with all data
print(f"Number of subjects with all data: {len(df)-missing_files_count}/{len(df)}")

# Convert the features and scores to arrays
X = np.array(features)
y = np.array(scores)
mask = ~np.isnan(y)
X = X[mask]
y = y[mask]


# Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

# Train the ridge regression model with cross-validation
reg = RidgeCV(alphas=[1e5,1e6,1e7,1e8], cv=5)
reg.fit(X_train, y_train)

y_pred=reg.predict(X_test)

# Calculate the correlation between predicted and observed values
corr = np.corrcoef(y_pred, y_test)[0, 1]
print(f"Correlation between predicted and observed values: {corr}")

# Permute the rows of the features matrix and predict on the permuted data
n_permutations = 100
null_corrs = []
for i in range(n_permutations):
    X_permuted = np.random.permutation(X_test)
    y_permuted_pred = reg.predict(X_permuted)
    null_corr = np.corrcoef(y_permuted_pred, y_test)[0, 1]
    null_corrs.append(null_corr)

# Calculate the p-value as the proportion of null correlations that are greater than the observed correlation
p_value = (np.sum(null_corrs > corr) + 1) / (n_permutations + 1)
print(f"P-value: {p_value}")

