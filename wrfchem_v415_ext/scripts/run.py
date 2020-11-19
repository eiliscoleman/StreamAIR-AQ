#!/usr/bin/env python3
import json
import subprocess as s


def main(config_path):
    config = {}
    with open(config_path) as f_config:
        config = json.load(f_config)

    for execute in config['execute']:
        s.run(config[execute]['program'],
              cwd=config[execute]['dir'],
              check=config[execute]['check'],
              shell=True)


if __name__ == '__main__':
    main('./config.json')
