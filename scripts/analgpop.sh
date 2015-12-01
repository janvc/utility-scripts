#!/bin/bash

rdgpop85 0 | sed -e '1,/Maximal/d' | sed -e '/-----/,$d' | sort -n

exit 0

