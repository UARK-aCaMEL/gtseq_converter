import pandas as pd

df = pd.read_csv("wtd_final_test/tables/wtd_test_sample_concordance_summary.csv")

# Define the threshold for high missing data percentage
threshold = 60.0

# Filter out rows with high missing data percentages
filtered_df = df[df["missing_pct"] <= threshold]

# Recalculate the concordant percentage
adjusted_concordant_pct = (
    filtered_df["concordant"].sum() / filtered_df["total"].sum()
) * 100

print(
    f"Adjusted Concordant Percentage (excluding high missing data): {adjusted_concordant_pct:.2f}%"
)

print(len(filtered_df))

print(1.96 * filtered_df["concordant_pct"].std())
