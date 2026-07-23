export USER=$(whoami)
export COBREMS_WORKER=/nfs/direct/packages/gpfs/gpfs1/osgusers/spotfinder
export ROOTSYS=//nfs/direct/packages/gpfs/gpfs1/osgusers/root_install
export PATH=$PATH:/usr/local/bin:$ROOTSYS/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib
export PYTHONPATH=$PYTHONPATH:$ROOTSYS/lib
export RABBITMQ_SERVER=nod30.phys.uconn.edu
export RABBITMQ_USER=spotfinder:FEA04D199530F7F723AE810
export RABBITMQ_VHOST=/spotfind
export REDIS_SERVER=nod30.phys.uconn.edu
export REDIS_PASSWORD="PCFK0yvDanTh00b6UHLZDNAW79527qmWIYjjZfwKvyg="

# start firewall services
#systemctl enable firewalld
#systemctl enable fail2ban
#systemctl restart firewalld
#systemctl restart fail2ban

# disable firmware GRO on em3 until I can figure out how to make it work
#ethtool -K em3 gro off
