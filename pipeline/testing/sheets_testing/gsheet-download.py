import json
import gspread
from oauth2client.service_account import ServiceAccountCredentials
import pandas
 
# use creds to create a client to interact with the Google Drive API
scope = ['https://spreadsheets.google.com/feeds']
creds = ServiceAccountCredentials.from_json_keyfile_name('Free Genes-153172d39301.json', scope)
client = gspread.authorize(creds)
 
# Find a workbook by name and open the first sheet
# Make sure you use the right name here.
sheet = client.open("Free Genes Project Rolling Submissions (Responses)").sheet1
 
# Extract and print all of the values
data = pandas.DataFrame(sheet.get_all_records())
print(data)

for index, row in data.iterrows():
    genbank = (row['Genbank file'])
    print (genbank)


