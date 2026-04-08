import datetime
import numpy as np
import pandas as pd


def parse_log_file(log_file: str):

    with open(log_file) as f:
        lines = f.readlines()
    lines = [line for line in lines if "cnt: " in line]

    logs = []
    for line in lines:
        log = {}
        token = line.split(" ")

        log["timestamp"] = datetime.datetime.strptime(token[0], "%Y-%m-%d,%H:%M:%S.%f")
        log["rig"] = token[1].rstrip(":")

        # for k, v in
        for k, v in zip(token[2:-1:2], token[3:-1:2]):
            k = k.rstrip(":")
            v = v.rstrip(";")
            try:
                log[k] = eval(v)
            except:
                log[k] = v

        logs.append(log)

    return pd.DataFrame(logs)
