#!/bin/bash
# script for cron to update docs on server
# 0 23 * * * cd /home/kelidas/workspace_git/partmod/docs; bash autobuild_docs.sh > /tmp/partmod.log

ACTION='\033[1;90m'
FINISHED='\033[1;96m'
READY='\033[1;92m'
NOCOLOR='\033[0m' # No Color
ERROR='\033[0;31m'

date
echo -e ${ACTION}Checking Git repo
echo -e =======================${NOCOLOR}
BRANCH=$(git rev-parse --abbrev-ref HEAD)
if [ "$BRANCH" != "master" ]
then
echo -e ${ERROR}Not on master. Aborting. ${NOCOLOR}
echo
exit 0
fi

git fetch
HEADHASH=$(git rev-parse HEAD)
UPSTREAMHASH=$(git rev-parse master@{upstream})

if [ "$HEADHASH" != "$UPSTREAMHASH" ]
then
    echo -e ${ERROR}Not up to date with origin. Aborting.${NOCOLOR}
    echo
    git pull

    cd ../../partmod-build
    cmake .
    cmake --build . --target docs
    echo -e ${FINISHED}Docs generated.${NOCOLOR}
    cmake --build . --target docs_publish
    echo -e ${FINISHED}Docs pushed.${NOCOLOR}
    cmake --build .
    echo -e ${FINISHED}DiscreteModel build.${NOCOLOR}

    build_win_dir="../partmod-build-win"
    if [ -d "$build_win_dir" ]
    then
        cd "$build_win_dir"
        cmake .
        cmake --build .
        echo -e ${FINISHED}DiscreteModel build for Win.${NOCOLOR}
    fi
else
    echo -e ${FINISHED}Current branch is up to date with origin/master.${NOCOLOR}
fi


