#!/bin/bash

password='rqliu1578'
sshpass -p $password scp -P 10190 rqliu@lssc-vpn:hsl/EFEDG/data/*.data ~/projects/EFEDG/data/.

