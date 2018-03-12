import pandas as pd
import json
import glob
from config import *

def query_db(query="*"):
    '''
    Pulls all of the current entries in the database and compiles them
    into a DataFrame. Creates compounding column names all starting with "data" to allow for
    nested keys to be reindexed later into their json files. Individual fragments are set as
    the index, i.e. a gene with two fragments would result in two different rows in the df.
    '''
    def crawl_dict(dictionary,name):
        '''
        Scans the dictionary to find the lowest level keys and their values, adds to a list of
        tuples and then returns the full list
        '''
        new_row = []
        for key in dictionary.keys():
            # The fragment sequences and locations columns are dictionaries
            # that have keys unique to the fragment name and so they need to
            # remain nested
            if key == "fragment_sequences" or key == "fragments":
                new_row = new_row + [(name+"."+key,dict(dictionary[key]))]
            elif type(dictionary[key]) == type(dict()):
                new_row += crawl_dict(dictionary[key],str(name+"."+key))
            elif type(dictionary[key]) == type(str()):
                new_row = new_row + [(name+"."+key,dictionary[key])]
            elif type(dictionary[key]) == type(list()) and dictionary[key] != [""]:
                new_row += crawl_dict(dictionary[key][-1],str(name+"."+key))
        return new_row

    # Iterates through all of the json files in the database and concatenates
    # all of the rows together to create the final dataframe
    full_db = pd.DataFrame()
    for file in glob.glob(BASE_PATH + "/data/*/BBF10K_{}.json".format(query)):
        with open(file,"r") as json_file:
            data = json.load(json_file)
        new_row = dict(crawl_dict(data,"data"))
        new_df = pd.DataFrame(new_row)
        full_db = pd.concat([full_db,new_df])
    return full_db

def write_db(df):
    '''
    Takes in a data frame that has columns named "data.[something].[else]"
    and it uses those names to index into the json file corresponding to
    the rows gene_id. Pulls each json file, changes the attributes listed
    in the dataframe and rewrites the file.
    '''
    # Creates a list of the column names so that they can  be included in the
    # name of each column
    keys = df.columns.values
    for index, row in df.iterrows():
        with open(BASE_PATH + "/data/{}/{}.json".format(row["data.gene_id"],row["data.gene_id"]),"r") as json_file:
            data = json.load(json_file)
        for key,value in zip(keys,row):
            # Creates a python command to access the desired location in the
            # nested dictionary
            k = key.split(".")
            location = "data"
            for item in k[1:]:
                location = location + '["{}"]'.format(item)
            # The fragment sequences and locations columns are dictionaries
            # that have keys unique to the fragment name and so those names
            # must be inserted
            if key == "data.sequence.fragment_sequences" or key == "data.location.fragments":
                location = location + '["{}"]'.format(index)
                exec('{} = "{}"'.format(location,value))
            else:
                exec('{} = "{}"'.format(location,value))
        with open(BASE_PATH + "/data/{}/{}.json".format(row["data.gene_id"],row["data.gene_id"]),"w+") as json_file:
            json.dump(data,json_file,indent=2)


full_db = query_db()
full_db.to_csv("./test_dict.csv")
write_db(full_db)
