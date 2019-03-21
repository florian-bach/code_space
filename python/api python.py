help=6
print("hello world")

# authenticate

import requests

url = "https://premium.cytobank.org/cytobank/api/v1/authenticate"

payload = "{ \"username\": \"emilsinclair\", \"password\": \"Laila1tov!123\" }"
headers = {
    'accept': "application/json",
    'charset': "UTF-8",
    'content-type': "application/json"
    }

response = requests.request("POST", url, data=payload, headers=headers)

print(response.text)

AUTH_TOKEN="eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJqdGkiOiJmM2UzMTE3Y2ExMGNlNzAzODdhY2I3M2E5NGNjNjhjYyIsImV4cCI6MTUzNDE4ODM2NCwidXNlcl9pZCI6NjYwMSwiYXVkIjoiY3l0b2JhbmtfYXBpX3YxX3VzZXJzIiwiaWF0IjoxNTM0MTU5NTY0LCJpc3MiOiJodHRwczovL3ByZW1pdW0uY3l0b2Jhbmsub3JnLyIsIm5iZiI6MTUzNDE1OTU2NCwic3ViIjoiY3l0b2JhbmtfYXBpX3YxIn0.g44mJMhYyGG5PCgjQaE7p02tUqtkTWPZm1RtbD6CPEg"

#create spade


url = "https://premium.cytobank.org/cytobank/api/v1/experiments/172509/advanced_analyses/spade/"

payload = "{ \"spade\": { \"name\":  \"APITEST\" } }"
headers = {
    'authorization': "Bearer eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJqdGkiOiJmM2UzMTE3Y2ExMGNlNzAzODdhY2I3M2E5NGNjNjhjYyIsImV4cCI6MTUzNDE4ODM2NCwidXNlcl9pZCI6NjYwMSwiYXVkIjoiY3l0b2JhbmtfYXBpX3YxX3VzZXJzIiwiaWF0IjoxNTM0MTU5NTY0LCJpc3MiOiJodHRwczovL3ByZW1pdW0uY3l0b2Jhbmsub3JnLyIsIm5iZiI6MTUzNDE1OTU2NCwic3ViIjoiY3l0b2JhbmtfYXBpX3YxIn0.g44mJMhYyGG5PCgjQaE7p02tUqtkTWPZm1RtbD6CPEg",
    'content-type': "application/json"
    }

response = requests.request("POST", url, data=payload, headers=headers)

print(response.text)


#copy spade analysis settings

url = "https://premium.cytobank.org/cytobank/api/v1/experiments/172509/advanced_analyses/spade/40547/copy_settings"

payload = ""
headers = {'authorization': 'Bearer eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJqdGkiOiJmM2UzMTE3Y2ExMGNlNzAzODdhY2I3M2E5NGNjNjhjYyIsImV4cCI6MTUzNDE4ODM2NCwidXNlcl9pZCI6NjYwMSwiYXVkIjoiY3l0b2JhbmtfYXBpX3YxX3VzZXJzIiwiaWF0IjoxNTM0MTU5NTY0LCJpc3MiOiJodHRwczovL3ByZW1pdW0uY3l0b2Jhbmsub3JnLyIsIm5iZiI6MTUzNDE1OTU2NCwic3ViIjoiY3l0b2JhbmtfYXBpX3YxIn0.g44mJMhYyGG5PCgjQaE7p02tUqtkTWPZm1RtbD6CPEg'}

response = requests.request("POST", url, data=payload, headers=headers)

print(response.text)


#rename spade analaysis



payload = "{ \"spade\": { \"name\": \"spade70\" } }"
headers = {
    'authorization': "Bearer eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJqdGkiOiJmM2UzMTE3Y2ExMGNlNzAzODdhY2I3M2E5NGNjNjhjYyIsImV4cCI6MTUzNDE4ODM2NCwidXNlcl9pZCI6NjYwMSwiYXVkIjoiY3l0b2JhbmtfYXBpX3YxX3VzZXJzIiwiaWF0IjoxNTM0MTU5NTY0LCJpc3MiOiJodHRwczovL3ByZW1pdW0uY3l0b2Jhbmsub3JnLyIsIm5iZiI6MTUzNDE1OTU2NCwic3ViIjoiY3l0b2JhbmtfYXBpX3YxIn0.g44mJMhYyGG5PCgjQaE7p02tUqtkTWPZm1RtbD6CPEg",
    'content-type': "application/json"
    }

response = requests.request("PUT", url, data=payload, headers=headers)

print(response.text)

#show spade details

url = "https://premium.cytobank.org/cytobank/api/v1/experiments/172509/advanced_analyses/spade/40559/"

querystring = {"include_settings":"1"}

headers = {'authorization': 'Bearer eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJqdGkiOiJmM2UzMTE3Y2ExMGNlNzAzODdhY2I3M2E5NGNjNjhjYyIsImV4cCI6MTUzNDE4ODM2NCwidXNlcl9pZCI6NjYwMSwiYXVkIjoiY3l0b2JhbmtfYXBpX3YxX3VzZXJzIiwiaWF0IjoxNTM0MTU5NTY0LCJpc3MiOiJodHRwczovL3ByZW1pdW0uY3l0b2Jhbmsub3JnLyIsIm5iZiI6MTUzNDE1OTU2NCwic3ViIjoiY3l0b2JhbmtfYXBpX3YxIn0.g44mJMhYyGG5PCgjQaE7p02tUqtkTWPZm1RtbD6CPEg'}

response = requests.request("GET", url, headers=headers, params=querystring)

print(response.text)


#update spade settings

url = "https://premium.cytobank.org/cytobank/api/v1/experiments/172509/advanced_analyses/spade/40562/"

payload = "{ \"spade\": { \"name\": \"APITEST\", \"compensationId\": -2, \"targetNumberOfNodes\": 20, \"population\": 0, \"clusteringChannels\": [2328,2329], \"downSampledEventsTarget\": { \"percent\": 100 }, \"foldChangeGroups":[{"name":"802","fcsFiles":[{"id":2842613,"name":"802+9_cd4_viSNE.fcs","baseline":false},{"id":2842614,"name":"802-1_cd4_viSNE.fcs","baseline":true},{"id":2842644,"name":"Specimen_001_802_007 C+8_cd4_viSNE.fcs","baseline":false}]},{"name":"806","fcsFiles":[{"id":2842615,"name":"806+10 C+10_cd4_viSNE.fcs","baseline":false},{"id":2842616,"name":"806c+09_cd4_viSNE.fcs","baseline":false},{"id":2842617,"name":"806c-01_cd4_viSNE.fcs","baseline":false},{"id":2842645,"name":"Specimen_001_806_002 C+8_cd4_viSNE.fcs","baseline":false}]},{"name":"808","fcsFiles":[{"id":2842618,"name":"808+9_cd4_viSNE.fcs","baseline":false},{"id":2842619,"name":"808-1_cd4_viSNE.fcs","baseline":true},{"id":2842646,"name":"Specimen_001_808_008 C+8_cd4_viSNE.fcs","baseline":false}]},{"name":"812","fcsFiles":[{"id":2842620,"name":"812+9_cd4_viSNE.fcs","baseline":false},{"id":2842621,"name":"812-1_cd4_viSNE.fcs","baseline":true},{"id":2842647,"name":"Specimen_001_812_009 C+8_cd4_viSNE.fcs","baseline":false}]},{"name":"814","fcsFiles":[{"id":2842622,"name":"814c-1_cd4_viSNE.fcs","baseline":true},{"id":2842648,"name":"Specimen_001_814_010 C+8_cd4_viSNE.fcs","baseline":false}]},{"name":"815","fcsFiles":[{"id":2842623,"name":"815+9_cd4_viSNE.fcs","baseline":false},{"id":2842624,"name":"815-1_cd4_viSNE.fcs","baseline":true},{"id":2842649,"name":"Specimen_001_815_011 C+8_cd4_viSNE.fcs","baseline":false}]},{"name":"818","fcsFiles":[{"id":2842625,"name":"818+9_cd4_viSNE.fcs","baseline":false},{"id":2842626,"name":"818-1_cd4_viSNE.fcs","baseline":true},{"id":2842650,"name":"Specimen_001_818_012 C+8_cd4_viSNE.fcs","baseline":false}]},{"name":"819","fcsFiles":[{"id":2842627,"name":"819+9_cd4_viSNE.fcs","baseline":false},{"id":2842628,"name":"819-1_cd4_viSNE.fcs","baseline":true},{"id":2842651,"name":"Specimen_001_819_013 C+8_cd4_viSNE.fcs","baseline":false}]},{"name":"820","fcsFiles":[{"id":2842629,"name":"820+9_cd4_viSNE.fcs","baseline":false},{"id":2842630,"name":"820-1_cd4_viSNE.fcs","baseline":true},{"id":2842652,"name":"Specimen_001_820_015 C+8_cd4_viSNE.fcs","baseline":false}]},{"name":"821","fcsFiles":[{"id":2842631,"name":"821+10 C+10_cd4_viSNE.fcs","baseline":false},{"id":2842632,"name":"821+9_cd4_viSNE.fcs","baseline":false},{"id":2842633,"name":"821-1_cd4_viSNE.fcs","baseline":true},{"id":2842653,"name":"Specimen_001_821_014 C+8_cd4_viSNE.fcs","baseline":false}]},{"name":"822","fcsFiles":[{"id":2842634,"name":"822+9_cd4_viSNE.fcs","baseline":false},{"id":2842635,"name":"822-1_cd4_viSNE.fcs","baseline":true},{"id":2842654,"name":"Specimen_001_822_016 C+8_cd4_viSNE.fcs","baseline":false}]},{"name":"823","fcsFiles":[{"id":2842636,"name":"823+10 C+10_cd4_viSNE.fcs","baseline":false},{"id":2842637,"name":"823+9_cd4_viSNE.fcs","baseline":false},{"id":2842638,"name":"823-1_cd4_viSNE.fcs","baseline":true},{"id":2842655,"name":"Specimen_001_823_017 C+8_cd4_viSNE.fcs","baseline":false}]},{"name":"824","fcsFiles":[{"id":2842639,"name":"824+10 C+10_cd4_viSNE.fcs","baseline":false},{"id":2842640,"name":"824+9_cd4_viSNE.fcs","baseline":false},{"id":2842641,"name":"824-1_cd4_viSNE.fcs","baseline":true},{"id":2842656,"name":"Specimen_001_824_018 C+8_cd4_viSNE.fcs","baseline":false}]},{"name":"825","fcsFiles":[{"id":2842642,"name":"825+9_cd4_viSNE.fcs","baseline":false},{"id":2842643,"name":"825-1_cd4_viSNE.fcs","baseline":true},{"id":2842657,"name":"Specimen_001_825_019 C+8_cd4_viSNE.fcs","baseline":false}]}]}
    'authorization': "Bearer eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJqdGkiOiJmM2UzMTE3Y2ExMGNlNzAzODdhY2I3M2E5NGNjNjhjYyIsImV4cCI6MTUzNDE4ODM2NCwidXNlcl9pZCI6NjYwMSwiYXVkIjoiY3l0b2JhbmtfYXBpX3YxX3VzZXJzIiwiaWF0IjoxNTM0MTU5NTY0LCJpc3MiOiJodHRwczovL3ByZW1pdW0uY3l0b2Jhbmsub3JnLyIsIm5iZiI6MTUzNDE1OTU2NCwic3ViIjoiY3l0b2JhbmtfYXBpX3YxIn0.g44mJMhYyGG5PCgjQaE7p02tUqtkTWPZm1RtbD6CPEg",
    'content-type': "application/json"
    }

response = requests.request("PUT", url, data=payload, headers=headers)

print(response.text)



#run spade



url = "https://premium.cytobank.org/cytobank/api/v1/experiments/172509/advanced_analyses/spade/40552/run"

payload = ""
headers = {'authorization': 'Bearer {{AUTH_TOKEN}}'}

response = requests.request("POST", url, data=payload, headers=headers)

print(response.text)
