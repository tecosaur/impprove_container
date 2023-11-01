FROM julia:1.9.3
LABEL maintainer="Timothy <timothy.chapman@telethonkids.org.au>"

RUN apt-get update \
    && apt-get -y install --no-install-recommends \
        gfortran gcc g++ \
    && rm -rf /var/lib/apt/lists/*

# We will put everything under the root directory '/'
# for compatibility with Apptainer.

COPY Project.toml Project.toml
COPY Manifest.toml Manifest.toml
COPY vendored-packages/IMPPROVE IMPPROVE
COPY vendored-packages/GeneVariantsMinimalData GeneVariantsMinimalData
RUN ln -s /data GeneVariantsMinimalData/data

ENV JULIA_DEPOT_PATH="/usr/local/julia/local/share/julia:/usr/local/julia/share/julia"
ENV DATATOOLKIT_STORE="/dtk_cache"
RUN mkdir -p /dtk_cache/store /dtk_cache/cache
RUN echo "inventory_version = 0\n" > /dtk_cache/Inventory.toml
RUN chmod -R a+r /dtk_cache

RUN julia --project=. -e 'using Pkg; \
    Pkg.instantiate(); \
    Pkg.develop(path="./IMPPROVE"); \
    Pkg.develop(path="./GeneVariantsMinimalData"); \
    Pkg.precompile()'

RUN julia --project="" -e 'using Pkg; \
    Pkg.activate(); \
    Pkg.add("PackageCompiler"); \
    using PackageCompiler; \
    Pkg.activate("."); \
    create_sysimage(collect(keys(Pkg.project().dependencies)), \
                    sysimage_path="sysimage.so", \
                    cpu_target = "generic;sandybridge,-xsaveopt,clone_all;haswell,-rdrnd,base(1)")'

# See <https://github.com/docker-library/julia/issues/79> for why this
# particular CPU target.

ADD https://github.com/brentp/slivar/releases/download/v0.3.0/slivar slivar
RUN chmod +x slivar
RUN cp slivar /usr/local/bin

COPY src/bcftools.jl bcftools.jl
RUN julia --project=/ bcftools.jl
RUN chmod +x bcftools
RUN cp bcftools /usr/local/bin

COPY src/apply.jl apply.jl
COPY src/apply_lib.jl apply_lib.jl
COPY src/render_lib.jl render_lib.jl
RUN _JL_DOCKER_PRERUN=1 julia -Jsysimage.so --project=. apply.jl

COPY README.txt README
COPY src/mini-readme.txt mini-readme.txt
COPY src/sample.vcf sample.vcf

COPY src/entrypoint.sh entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
