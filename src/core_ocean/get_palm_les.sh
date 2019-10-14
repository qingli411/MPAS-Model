#!/bin/bash

## LES Tag for build
LES_TAG=25c5e73
## Subdirectory in LES repo to use
LES_SUBDIR=trunk/SOURCE

## Available protocols for acquiring LES source code
LES_GIT_HTTP_ADDRESS=https://github.com/qingli411/palm_les_lanl.git
LES_GIT_SSH_ADDRESS=git@github.com:qingli411/palm_les_lanl.git

GIT=`which git`
PROTOCOL=""

# LES exists. Check to see if it is the correct version.
# Otherwise, flush the directory to ensure it's updated.
if [ -d palm_les ]; then

	if [ -d .les_all/.git ]; then
		cd .les_all
		CURR_TAG=$(git rev-parse --short HEAD)
		cd ../
		if [ "${CURR_TAG}" == "${LES_TAG}" ]; then
			echo "PALM LES version is current. Skip update"
		else
			unlink palm_les
			rm -rf .les_all
		fi
	else
		unlink palm_les
		rm -rf .les_all
	fi
fi


# LES Doesn't exist, need to acquire souce code
# If might have been flushed from the above if, in the case where it was svn or wget that acquired the source.
if [ ! -d palm_les ]; then
	if [ -d .les_all ]; then
		rm -rf .les_all
	fi

		echo " ** Using git to acquire LES source. ** "
		PROTOCOL="git ssh"
		git clone ${LES_GIT_SSH_ADDRESS} .les_all &> /dev/null
		if [ -d .les_all ]; then
			cd .les_all
			git checkout ${LES_TAG} &> /dev/null
			cd ../
			ln -sf .les_all/${LES_SUBDIR} palm_les
		else
			git clone ${LES_GIT_HTTP_ADDRESS} .les_all &> /dev/null
			PROTOCOL="git http"
			if [ -d .les_all ]; then
				cd .les_all
				git checkout ${LES_TAG} &> /dev/null
				cd ../
				ln -sf .les_all/${LES_SUBDIR} palm_les
			fi
		fi
fi

if [ ! -d palm_les ]; then
	echo " ****************************************************** "
	echo " ERROR: Build failed to acquire LES source."
	echo ""
	echo " Please ensure your proxy information is setup properly for"
	echo " the protocol you use to acquire LES."
	echo ""
	echo " The automated script attempted to use: ${PROTOCOL}"
	echo ""
	if [ "${PROTOCOL}" == "git http" ]; then
		echo " This protocol requires setting up the http.proxy git config option."
	elif [ "${PROTOCOL}" == "git ssh" ]; then
		echo " This protocol requires having ssh-keys setup, and ssh access to git@github.com."
		echo " Please use 'ssh -vT git@github.com' to debug issues with ssh keys."
	fi
	echo ""
	echo " ****************************************************** "
fi
