#!/bin/bash

if [[ $# -ne 1 ]]; then
    echo "usage: run-web-viewer.sh <port>"
    exit 1
fi

uwsgi --socket 0.0.0.0:$1 --plugin python --protocol=http -w web_viewer --callable app
