#!/usr/bin/env python
# coding: utf-8

# # Pymaceuticals Inc.
# ---
# 
# ### Analysis
# 
# * From the bar graph analysis, it can be said that maximum number of mice are tested for drugs 'Capomulin' and 'Ramicane'.
# * The tests performed on mice has 51% of male mice and 49% of female mice.
# * It is found that the final tumor volume of mice treated with 'Capomulin' and 'Ramicane' drugs was very less compared to final tumor volume of mice tested with 'Infubinol' and 'Ceftamin' drugs.(almost 20 mm3 less)
# * After treating mouse l509 with 'Capomulin' drug, there is a slight increase in tumor volume initially but after few timepoints, decrease in the tumour volume can be observed.
# * The correlation between mouse weight and the average tumor volume for the Capomulin regimen drug is 0.84. There is a positive corelation between mouse weight and the average tumor volume in this case.

# In[2]:


# Dependencies and Setup
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats as st

# Study data files
mouse_metadata_path = "data/Mouse_metadata.csv"
study_results_path = "data/Study_results.csv"

# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)

# Combine the data into a single dataset
combined_data = pd.merge(study_results,mouse_metadata, how = "left", on = "Mouse ID")
# Display the data table for preview
combined_data


# In[3]:


# Checking the number of mice.
total_mice = len(combined_data["Mouse ID"].unique())
total_mice


# In[4]:


# Getting the duplicate mice by ID number that shows up for Mouse ID and Timepoint. 
duplicate_ids = combined_data[combined_data.duplicated()]["Mouse ID"].unique()
duplicate_ids


# In[5]:


# Optional: Get all the data for the duplicate mouse ID. 
combined_data[combined_data ["Mouse ID"] == "g989"]


# In[6]:


# Create a clean DataFrame by dropping the duplicate mouse by its ID.
clean_data_df = combined_data[combined_data["Mouse ID"].isin(duplicate_ids) == False]
clean_data_df


# In[7]:


# Checking the number of mice in the clean DataFrame.
len(clean_data_df["Mouse ID"].unique())


# ## Summary Statistics

# In[8]:


# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen

# Use groupby and summary statistical methods to calculate the following properties of each drug regimen: 
# mean, median, variance, standard deviation, and SEM of the tumor volume. 
# Assemble the resulting series into a single summary DataFrame.
mean = clean_data_df.groupby("Drug Regimen").mean()["Tumor Volume (mm3)"]
median = clean_data_df.groupby("Drug Regimen").median()["Tumor Volume (mm3)"]
std_dev = clean_data_df.groupby("Drug Regimen").std()["Tumor Volume (mm3)"]
variance = clean_data_df.groupby("Drug Regimen").var()["Tumor Volume (mm3)"]
error = clean_data_df.groupby("Drug Regimen").sem()["Tumor Volume (mm3)"]

summary_df = pd.DataFrame(
    {"Mean Tumor Volume":mean,
    "Median Tumor Volume":median,
    "Tumor Volume Variance" : variance,
    "Tumor Volume Std. Dev.": std_dev,
    "Tumor Volume Std. Err.": error}
)
summary_df


# In[9]:


# Generate a summary statistics table of mean, median, variance, standard deviation, 
# and SEM of the tumor volume for each regimen

# Using the aggregation method, produce the same summary statistics in a single line.
summary = clean_data_df.groupby("Drug Regimen").agg({"Tumor Volume (mm3)":["mean","median","var","std","sem"]})
summary


# ## Bar and Pie Charts

# In[10]:


# Generate a bar plot showing the total number of timepoints for all mice tested for each drug regimen using Pandas.

mice_tested_counts = clean_data_df["Drug Regimen"].value_counts()
mice_tested_counts.plot(kind = "bar", figsize = (6,4))
plt.xlabel("Drug Regimen")
plt.ylabel("Number of Mice Tested")


# In[11]:


# Generate a bar plot showing the total number of timepoints for all mice tested for each drug regimen using pyplot.
mice_tested_counts = clean_data_df["Drug Regimen"].value_counts()
plt.bar(mice_tested_counts.index.values, mice_tested_counts.values)
plt.xticks(rotation = 90)
plt.xlabel("Drug Regimen")
plt.ylabel("Number of Mice Tested")
plt.show()


# In[12]:


# Generate a pie plot showing the distribution of female versus male mice using Pandas
mice_sex_counts = clean_data_df["Sex"].value_counts()
mice_sex_counts.plot(kind = "pie", autopct="%1.1f%%")


# In[13]:


# Generate a pie plot showing the distribution of female versus male mice using pyplot
mice_sex_counts = clean_data_df["Sex"].value_counts()
plt.pie(mice_sex_counts.values, autopct="%1.1f%%" , labels = mice_sex_counts.index.values)
plt.ylabel("Sex")


# ## Quartiles, Outliers and Boxplots

# In[14]:


# Calculate the final tumor volume of each mouse across four of the treatment regimens:  
# Capomulin, Ramicane, Infubinol, and Ceftamin

# Start by getting the last (greatest) timepoint for each mouse

max_tumor = clean_data_df.groupby("Mouse ID")["Timepoint"].max()
max_tumor = max_tumor.reset_index()

# Merge this group df with the original DataFrame to get the tumor volume at the last timepoint
merged_max_tumor = max_tumor.merge(clean_data_df, on = ["Mouse ID","Timepoint"], how = "left")
merged_max_tumor


# In[15]:


# Put treatments into a list for for loop (and later for plot labels)
treatments = ["Capomulin", "Ramicane","Infubinol","Ceftamin"]

# Create empty list to fill with tumor vol data (for plotting)
tumor_volume  = []

# Calculate the IQR and quantitatively determine if there are any potential outliers. 
for x in treatments:
    
    # Locate the rows which contain mice on each drug and get the tumor volumes
    final_tumor_volume = merged_max_tumor.loc[merged_max_tumor["Drug Regimen"] == x, "Tumor Volume (mm3)"]
    
    # add subset 
    
    tumor_volume.append(final_tumor_volume)
    # Determine outliers using upper and lower bounds
    quartiles = final_tumor_volume.quantile([0.25,0.5,0.75])
    lowerq = quartiles[0.25]
    upperq = quartiles[0.75]
    iqr =  upperq - lowerq
    lower_bound = lowerq - (1.5 * iqr)
    upper_bound = upperq + (1.5 * iqr)
    outliers = final_tumor_volume.loc[(final_tumor_volume < lower_bound) | (final_tumor_volume > upper_bound)]
    
    print(f"{x}'s potential outliers {outliers}")


# In[16]:


# Generate a box plot that shows the distrubution of the tumor volume for each treatment group.
format_plot = dict(markerfacecolor = 'red', markersize = 8)
plt.boxplot(tumor_volume, labels = treatments, flierprops = format_plot)
plt.ylabel("Final Tumor Volume (mm3)")


# ## Line and Scatter Plots

# In[17]:


# Generate a line plot of tumor volume vs. time point for a mouse treated with Capomulin

Capomulin_data = clean_data_df[clean_data_df["Drug Regimen"] == "Capomulin"]
mouse_data = clean_data_df[clean_data_df["Mouse ID"] == "l509"]
plt.plot(mouse_data["Timepoint"], mouse_data["Tumor Volume (mm3)"])
plt.ylabel("Tumor Volume (mm3)")
plt.xlabel("Timepoint (Days)")
plt.title("Capomulin Treatment of Mouse l509")


# In[159]:


# Generate a scatter plot of average tumor volume vs. mouse weight for the Capomulin regimen
Capomulin_data = clean_data_df[clean_data_df["Drug Regimen"] == "Capomulin"]
avg_data = Capomulin_data.groupby("Mouse ID").mean()
plt.scatter(avg_data["Weight (g)"], avg_data["Tumor Volume (mm3)"])
plt.xlabel("Weight (g)")
plt.ylabel("Average Tumor Volume (mm3)")


# ## Correlation and Regression

# In[167]:


# Calculate the correlation coefficient and linear regression model 
# for mouse weight and average tumor volume for the Capomulin regimen

correlation = st.pearsonr(avg_data["Weight (g)"], avg_data[("Tumor Volume (mm3)")])
print(f"The correlation between mouse weight and the average tumor volume is {round(correlation[0],2)}")

values = st.linregress(avg_data["Weight (g)"], avg_data[("Tumor Volume (mm3)")])
slope = values[0]
intercept = values[1]
y_values = avg_data["Weight (g)"] * slope + intercept
plt.scatter(avg_data["Weight (g)"], avg_data[("Tumor Volume (mm3)")])
plt.plot(avg_data["Weight (g)"], y_values, color ="red" )
plt.xlabel("Weight (g)")
plt.ylabel("Average Tumor Volume (mm3)")


# In[ ]:




