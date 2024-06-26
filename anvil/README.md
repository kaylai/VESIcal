README for using Anvil server to run a remote instance of thermoengine and VESIcal on Digital Ocean droplet
===========================================================================================================
If you try this guide and run into any issues, please let me know! Email me. kayla.iacovino@nasa.gov. Cheers!

Steps overview:
1. Install docker image of thermoengine on DigitalOcean droplet
2. Install dependencies on DO droplet
3. Create a new screen and run the docker container
4. Create a script to relaunch the anvil server python script if it's not already running
5. Create a cron job to run your bash script every minute

# 1. Install docker image of thermoengine on DigitalOcean droplet
1. Create new droplet on DigitalOcean: Ubuntu 18.04
	1. See https://www.digitalocean.com/docs/droplets/how-to/create/
2. SSH into droplet (or use console accessed on DigitalOcean portal)
	1. On your mac, open terminal
	2. Type ssh root@[ip-address-of-droplet] (without brackets)
	3. Enter your DO droplet password
3. Install dependencies on DO droplet (see https://www.digitalocean.com/community/tutorials/how-to-install-and-use-docker-on-ubuntu-18-04)
```
sudo apt update
sudo apt install apt-transport-https ca-certificates curl software-properties-common
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu bionic stable"
sudo apt update
apt-cache policy docker-ce
```
4. Install Docker on DO droplet:
```
sudo apt install docker-ce
sudo systemctl status docker
```
5. Allow docker commands without sudo:
```
sudo usermod -aG docker ${USER}
su - ${USER}
id -nG
```
6. Install the thermoengine docker image on the DO droplet:
```
docker pull registry.gitlab.com/enki-portal/thermoengine:latest
```

# 2. Install dependencies on DO droplet
1. Install anaconda3 on the DO droplet:
	1. ```wget "https://repo.anaconda.com/archive/Anaconda3-2020.02-Linux-x86_64.sh"```
	2. ```bash Anaconda3-2020.02-Linux-x86_64.sh```
2. Install screen and cron:
	1. ```sudo apt-get install cron```
		1. It'll probably ask you what editor you want to use. We recommend nano, as does cron
	2. ```sudo apt-get install screen```

# 3. Create a new screen and run the docker container
1. Create a new screen and give it a name:
	1. ```screen -S <name-of-your-screen>```
	2. Note: to detach from this screen at any time, press ctrl+a then ctrl+d
	3. Note: to reattach to this screen at any time, type:
		1. You can see names of screens with `screen -ls`
		2. ```screen -r <name-of-your-screen>```
2. In your screen, run the docker image and give it a name. Here we use the name 'anvilcontainer':
	1. ```docker run --name anvilcontainer --user root -it --rm registry.gitlab.com/enki-portal/thermoengine:latest bash```
3. In the docker container, clone the VESIcal repository as:
	1. ```git clone https://github.com/kaylai/VESIcal.git```
4. In the docker container, install the anvil server library and VESIcal:
	1. ```pip install anvil-app-server```
	2. ```pip install VESIcal```
5. Detach from your screen:
	1. ctrl+a then ctrl+d

# 4. Create a script to relaunch the anvil server python script if it's not already running
1. Use nano or your desired editor to create and start editing a file, let's call it launch_anvil_server.sh
	1. nano launch_anvil_server.sh
	2. Note: to save your file in nano use ctrl+o and then enter to confirm file name
	3. Note: to exit nano use ctrl+x
2. In that file put the following (here we use the container name anvilcontainer. You must use the name of the container that you created above):
```
#!/bin/bash
truncate -s 0 cronfile.log
echo $(date)
if pgrep -f "python anvil_server.py" &>/dev/null; then
        echo "it is running"
else
        echo "starting anvil server anew"
        docker exec -i -w /app/VESIcal anvilcontainer python anvil_server.py
fi
```

4. Note: "/app/VESIcal" is the path to the folder VESIcal cloned from the github repo (that lives inside the docker container). In the case of the thermoengine docker, you live in the /app folder, so the path is /app/VESIcal. You'll be able to see the /app folder if you are inside of the running docker container that you created, which should be named anvilcontainer. Do `docker ps` to see containers. Names are shown on the far right.

# 5. Create a cron job to run your bash script every minute
1. Open up your crontab for editing with:
```crontab -e```
2. In your crontab add this line:
```* * * * * <path-to-your_bash_script> >> <path-to-a-place-you-like>/cronfile.log 2>&1```
3. Save and exit crontab (nano ctrl+o, enter, ctrl+x; vim esc wq)
4. Tip: If you don't know where your bash script is, when not in the crontab, go to the directory with your bash script in it and type pwd. That's the path to that location.
