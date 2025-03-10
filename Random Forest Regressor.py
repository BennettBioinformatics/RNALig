import pandas as pd
from sklearn.ensemble import RandomForestRegressor 
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
import joblib
from sklearn.metrics import r2_score, mean_squared_error

# Load the uploaded training data
training_file_path = '.csv'  # Path to the uploaded training data
training_data = pd.read_csv(training_file_path)

# Define the features and target
features = [
    "nucleotide_composition_A", "nucleotide_composition_U",
    "nucleotide_composition_G", "nucleotide_composition_C",
    "gc_content", "minimum_free_energy", "Molecular Weight (g/mol)",
    "Hydrogen Bond Donors", "Hydrogen Bond Acceptors", "LogP", "TPSA (Å²)",
    "Atom Count", "Chiral Atom Count", "num_electrostatic_contacts",
    "avg_electrostatic_distance", "num_vdw_contacts", "avg_vdw_distance"
]
target = "Binding Affinity (kcal/mol)"

# Split the data into features (X) and target (y)
X_train = training_data[features]
y_train = training_data[target]

# Step 2: Preprocess the data
imputer = SimpleImputer(strategy='mean')  # For missing values imputation
scaler = StandardScaler()  # For feature scaling
target_scaler = StandardScaler()  # For target scaling

# Handle missing values and scale features
X_train_imputed = imputer.fit_transform(X_train)  # Apply imputation on training features
X_train_scaled = scaler.fit_transform(X_train_imputed)  # Scale the features

# Scale target
y_train_scaled = target_scaler.fit_transform(y_train.values.reshape(-1, 1)).ravel()  # Reshaping target for scaling

# Step 3: Initialize Random Forest Regressor and train it
rf_model = RandomForestRegressor(random_state=42)

# Hyperparameter tuning with GridSearchCV (optional)
param_grid_rf = {
    'n_estimators': [100, 200, 500],
    'max_depth': [5, 10, 20, None],
    'min_samples_split': [5, 10, 20],
    'min_samples_leaf': [2, 4, 6],
    'max_features': ['sqrt', 'log2', None]
}

grid_search_rf = GridSearchCV(estimator=rf_model, param_grid=param_grid_rf, cv=5, scoring='r2', n_jobs=-1)
grid_search_rf.fit(X_train_scaled, y_train_scaled)

# Get the best model from GridSearchCV
best_rf_model = grid_search_rf.best_estimator_
print(f"Best Parameters (RF): {grid_search_rf.best_params_}")

# Step 4: Save the trained model and preprocessors
joblib.dump(imputer, 'imputer.pkl')  # Save the imputer (for missing value imputation)
joblib.dump(scaler, 'scaler.pkl')  # Save the feature scaler
joblib.dump(target_scaler, 'target_scaler.pkl')  # Save the target scaler
joblib.dump(best_rf_model, 'best_rf_model.pkl')  # Save the trained Random Forest model

print("Model and preprocessing objects saved.")

# Step 5: Predict Binding Affinity for New PDB

def predict_binding_affinity(new_pdb_file):
    # Load the new PDB data with features (excluding the target column)
    new_data = pd.read_csv(new_pdb_file)
    
    # Ensure the new data contains the required features
    if not all(feature in new_data.columns for feature in features):
        raise ValueError("The new data does not contain the required features.")

    # Extract the features for prediction
    X_new = new_data[features]

    # Preprocess the new data (imputation and scaling)
    X_new_imputed = imputer.transform(X_new)
    X_new_scaled = scaler.transform(X_new_imputed)

    # Predict the binding affinity using the trained model
    y_new_pred_scaled = best_rf_model.predict(X_new_scaled)

    # Inverse transform the predictions to the original scale
    y_new_pred = target_scaler.inverse_transform(y_new_pred_scaled.reshape(-1, 1)).flatten()

    # Save the predictions to a new CSV file
    new_pdb_results = pd.DataFrame({
        "PDB ID": new_data["PDB ID"],
        "Predicted Binding Affinity": y_new_pred
    })
    
    # Save the results to CSV
    output_filename = f"Predicted_Binding_Affinity_{new_pdb_file.split('/')[-1]}"
    new_pdb_results.to_csv(output_filename, index=False)
    print(f"Predictions saved to '{output_filename}'")

# Step 6: User Input to Predict Binding Affinity for New PDB

# Get the path of the new PDB features CSV file from the user
new_pdb_file = input("Enter the path of the new PDB features CSV file: ")
predict_binding_affinity(new_pdb_file)
