#!/bin/bash

set -u

TARGET_HEADER="${1}"

REV_TEMPLATE="@COMMIT@"

#don't fail if no .git
git branch &> /dev/null || exit 0

REV="$(git describe --tags)"
NEWLINE='#define COMMIT "'${REV}'"'
echo $NEWLINE

sed -i "s/^#define.*COMMIT.*$/${NEWLINE}/g"  ${TARGET_HEADER}


