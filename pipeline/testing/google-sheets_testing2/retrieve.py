import requests
response = requests.get('https://docs.google.com/spreadsheet/ccc?key=1rR4RQ1GLN3eMkOHVRc3jXZ6ZR2FfbUGU9JUXkm8OFHk&output=csv')
assert response.status_code == 200, 'Wrong status code'
print (response.content)

https://docs.google.com/spreadsheets/d/1rR4RQ1GLN3eMkOHVRc3jXZ6ZR2FfbUGU9JUXkm8OFHk/edit?usp=sharing

