#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process ManhattanPlot {
    input:

    output:
        path(manhattan)
    
    script:
        """
        set -eo pipefail

        
        """
}