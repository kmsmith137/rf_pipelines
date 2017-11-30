#!/bin/bash

if [[ $# -ne 2 ]]; then
    echo "Usage: run-web-viewer.sh <root_dir> <port>"
    echo "On frb1, you probably want: run-web-viewer.sh /data2/web_viewer 5000"
    echo "On cf0g9, you probably want: run-web-viewer.sh /frb-archiver-1/web_viewer 5003"
    exit 1
fi

export WEB_VIEWER_ROOT=$1
uwsgi --socket 0.0.0.0:$2 --processes 8 --plugin python --protocol=http -w web_viewer --callable app
