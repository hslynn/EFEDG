#!/bin/bash

select isvpn in "VPN" "No VPN"; do
    case $isvpn in
        "VPN" ) ipname=lssc-vpn; break;;
        "No VPN" ) ipname=lssc; break;;
    esac
done

password='rqliu1578'
sshpass -p $password scp -P 10190 $* rqliu@$ipname:hsl/
