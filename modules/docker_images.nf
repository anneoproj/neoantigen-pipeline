process ENSURE_DOCKER_IMAGES {

    tag "docker-images"
    cache false

    output:
    path("docker_images.ready")

    script:
    """
    set -euo pipefail

    docker build --platform linux/amd64 -t netchop_image:latest ${projectDir}/docker/netchop

    docker build --platform linux/arm64/v8 -t netmhcpan_image:latest ${projectDir}/docker/netmhcpan

    docker build --platform linux/amd64 -t mhcflurry_image:latest ${projectDir}/docker/mhcflurry

    docker build --platform linux/amd64 -t mixmhcpred_image:latest ${projectDir}/docker/mixmhcpred

    docker image inspect netchop_image:latest --format '{{.Architecture}}' | grep -qx amd64
    docker image inspect netmhcpan_image:latest --format '{{.Architecture}}' | grep -qx arm64
    docker image inspect mhcflurry_image:latest --format '{{.Architecture}}' | grep -qx amd64
    docker image inspect mixmhcpred_image:latest --format '{{.Architecture}}' | grep -qx amd64

    touch docker_images.ready
    """
}

process ENSURE_METRIC_DOCKER_IMAGES {

    tag "metric-docker-images"
    cache false

    output:
    path("metric_docker_images.ready")

    script:
    """
    set -euo pipefail

    # NOTE: legacy workflow expects a separate channel-ready signal for
    # metrics branch. Kept as a dedicated target for compatibility.
    test -f docker_images.ready || (
        docker build --platform linux/amd64 -t netchop_image:latest ${projectDir}/docker/netchop
        docker build --platform linux/arm64/v8 -t netmhcpan_image:latest ${projectDir}/docker/netmhcpan
        docker build --platform linux/amd64 -t mhcflurry_image:latest ${projectDir}/docker/mhcflurry
        docker build --platform linux/amd64 -t mixmhcpred_image:latest ${projectDir}/docker/mixmhcpred
        docker image inspect netchop_image:latest --format '{{.Architecture}}' | grep -qx amd64
        docker image inspect netmhcpan_image:latest --format '{{.Architecture}}' | grep -qx arm64
        docker image inspect mhcflurry_image:latest --format '{{.Architecture}}' | grep -qx amd64
        docker image inspect mixmhcpred_image:latest --format '{{.Architecture}}' | grep -qx amd64
        touch docker_images.ready
    )
    touch metric_docker_images.ready
    """
}

process ENSURE_NEOPEP_DOCKER_IMAGES {

    tag "neopep-docker-images"
    cache false

    output:
    path("neopep_docker_images.ready")

    script:
    """
    set -euo pipefail

    # NOTE: placeholder compatibility step for neopep workflows.
    test -f docker_images.ready || (
        docker build --platform linux/amd64 -t netchop_image:latest ${projectDir}/docker/netchop
        docker build --platform linux/arm64/v8 -t netmhcpan_image:latest ${projectDir}/docker/netmhcpan
        docker build --platform linux/amd64 -t mhcflurry_image:latest ${projectDir}/docker/mhcflurry
        docker build --platform linux/amd64 -t mixmhcpred_image:latest ${projectDir}/docker/mixmhcpred
        docker image inspect netchop_image:latest --format '{{.Architecture}}' | grep -qx amd64
        docker image inspect netmhcpan_image:latest --format '{{.Architecture}}' | grep -qx arm64
        docker image inspect mhcflurry_image:latest --format '{{.Architecture}}' | grep -qx amd64
        docker image inspect mixmhcpred_image:latest --format '{{.Architecture}}' | grep -qx amd64
        touch docker_images.ready
    )
    touch neopep_docker_images.ready
    """
}
