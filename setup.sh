export USER=$(whoami)
export COBREMS_WORKER=/home/$USER/spotfinder
export PATH=$PATH:/usr/local/bin:/cvmfs/oasis.opensciencegrid.org/gluex/root-6.22.06/x86_64/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cvmfs/oasis.opensciencegrid.org/gluex/root-6.22.06/x86_64/lib
export PYTHONPATH=$PYTHONPATH:/cvmfs/oasis.opensciencegrid.org/gluex/root-6.22.06/x86_64/lib
export RABBITMQ_SERVER=cn410.storrs.hpc.uconn.edu
export RABBITMQ_USER=spotfinder:FEA098B2D91AC7F723AE810
export RABBITMQ_VHOST=spotfind

# start firewall services
#systemctl enable firewalld
#systemctl enable fail2ban
#systemctl restart firewalld
#systemctl restart fail2ban

# disable firmware GRO on em3 until I can figure out how to make it work
#ethtool -K em3 gro off
