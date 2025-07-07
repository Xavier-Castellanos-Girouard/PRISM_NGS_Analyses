#!/usr/bin/env bash
for f in `curl -k "https://genomique.iric.ca/BamList?key=2409-75d8a78bf3b74bceeefde0b34bf093c6&projectID=1374"`
do
  curl -k -o "#1" -C - --create-dirs -O $f
done
