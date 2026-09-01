#!/bin/bash
sudo docker pull ghcr.io/rprebet/albatross-2026:latest
sudo docker run --rm -it -p 8881:8881 --mount type=bind,src=".",dst=/home/vscode ghcr.io/rprebet/albatross-2026:latest
