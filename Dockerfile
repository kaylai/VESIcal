FROM registry.gitlab.com/enki-portal/thermoengine:master

# Make sure the contents of our repo are in ${HOME}
COPY . ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
RUN pip install --no-cache-dir appmode
RUN pip install --no-cache-dir VESIcal
RUN jupyter nbextension enable --py --sys-prefix appmode
RUN jupyter serverextension enable --py --sys-prefix appmode
RUN mkdir ThermoEngine # make this directory
USER ${NB_USER}

# create a user called 1000 (binder cannot run as root)
# ARG NB_USER=jovyan
# ARG NB_UID=1000
# ENV USER ${NB_USER}
# ENV NB_UID ${NB_UID}
# ENV HOME /home/${NB_USER}

# RUN adduser --disabled-password \
#     --gecos "Default user" \
#     --uid ${NB_UID} \
#     ${NB_USER}