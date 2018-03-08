import json
import gspread
from oauth2client.service_account import ServiceAccountCredentials
import pandas as pd
import sys
import numpy as np


# use creds to create a client to interact with the Google Drive API
scope = ['https://spreadsheets.google.com/feeds']
creds = ServiceAccountCredentials.from_json_keyfile_name('Free Genes-153172d39301.json', scope)
client = gspread.authorize(creds)

# Find a workbook by name and open the first sheet
sheet = client.open("Free Genes Project Rolling Submissions (Responses)").sheet1

# Extract all of the values
data = pd.DataFrame(sheet.get_all_records())
data = data.fillna('')
previous = pd.read_csv('Free Genes Project Rolling Submissions (Responses) - Form Responses 1.csv')
previous = previous.fillna('')
# Extract unique data
united_data = pd.concat([data, previous])
print(united_data)

united_data_grouped = united_data.groupby(list(united_data.columns))
uniq_data_idx = [x[0] for x in united_data_grouped.indices.values() if len(x) == 1]
uniq_data = united_data.iloc[uniq_data_idx]
print(uniq_data)

#print (uniq_data)

# Next id in the database
def NextID():
    data = (sorted(glob.glob("./data/*")))
    string = (data[-1])
    number = int(string[-6:]) + 1
    string_number = str(number)
    id_number = (string_number.zfill(6))
    full_id = "BBF10K_" + id_number
    return full_id
