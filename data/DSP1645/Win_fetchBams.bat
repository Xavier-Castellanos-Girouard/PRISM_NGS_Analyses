
curl -k "https://genomique.iric.ca/BamList?key=2409-75d8a78bf3b74bceeefde0b34bf093c6&projectID=1645" > 1645.bamlist
for /F "tokens=*" %%A in (1645.bamlist) do curl -k -o "#1" -C - --create-dirs -O %%A 
