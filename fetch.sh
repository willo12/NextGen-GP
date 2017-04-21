#!/usr/bin/env bash

name="$1"

ssh katana "cd ${kpath};tar cvzf export.tar.gz $1/{elite*,report*,config.ini}";
scp katana:${kpath}export.tar.gz ./;
tar xvzf export.tar.gz
