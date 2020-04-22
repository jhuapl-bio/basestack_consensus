#!/bin/bash

sed '$!N;s/\n/\t/' "$@"
