
curl -k "https://genomique.iric.ca/FastQList?key=2409-75d8a78bf3b74bceeefde0b34bf093c6&projectID=1472" > 1472.filelist
for /F "tokens=*" %%A in (1472.filelist) do curl -k -o "#1" -C - --create-dirs -O %%A
