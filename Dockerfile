FROM 812206152185.dkr.ecr.us-west-2.amazonaws.com/13502_wf_init_coexpression_wf:0.7.6-c0086f

RUN apt-get update -y && \
    apt-get install -y \
        libgdal-dev \
        libmagick++-dev

# Copy files for .renvignore to work
RUN rm /root/*.Rproj
COPY neighborhood_wf.Rproj /root/neighborhood_wf.Rproj
COPY .renvignore /root/.renvignore
COPY download_git.R /root/download_git.R
RUN Rscript /root/download_git.R
COPY renv.lock /root/renv.lock
RUN R -e "renv::restore(prompt = FALSE)"

# STOP HERE:
# The following lines are needed to ensure your build environement works
# correctly with latch.
RUN python3 -m pip install --upgrade latch
COPY wf /root/wf
ARG tag
ENV FLYTE_INTERNAL_IMAGE $tag
WORKDIR /root
