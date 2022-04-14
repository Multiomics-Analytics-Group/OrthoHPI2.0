import React, { FC } from "react";
import { BsInfoCircle } from "react-icons/bs";

import Panel from "./Panel";

const DescriptionPanel: FC = () => {
  return (
    <Panel
      initiallyDeployed
      title={
        <>
          <BsInfoCircle className="text-muted" /> Description
        </>
      }
    >
      <p>
        This graph represents <i>Predicted Human-parasite Protein-Protein Interactions</i> for 20 parasitic species. Each{" "}
        <i>node</i> represents a protein, and each edge a predicted interaction.
      </p>
      <p>
        The interactions are predicted based on homology extracted from{" "}
        <a target="_blank" rel="noreferrer" href="http://eggnog5.embl.de/">
          EggNOG 5.0
        </a> and intra-species from{" "}
        <a target="_blank" rel="noreferrer" href="https://string-db.org">
          STRING 11.5
        </a>
        , and then filtered using tissue information from{" "}
        <a target="_blank" rel="noreferrer" href="https://tissues.jensenlab.org">
          TISSUES
        </a>
        , cellular compartments from{" "}
        <a target="_blank" rel="noreferrer" href="https://compartments.jensenlab.org">
          COMPARTMENTS
        </a>
        , and secretome predictions from{" "}
        <a target="_blank" rel="noreferrer" href="https://wikipedia.org">
          DeepPredictions
        </a>
        .
      </p>
      <p>
        This web application has been developed by the Multi-omics Graph Analytics Group
        , using{" "}
        <a target="_blank" rel="noreferrer" href="https://reactjs.org/">
          react
        </a>{" "}
        and{" "}
        <a target="_blank" rel="noreferrer" href="https://www.sigmajs.org">
          sigma.js
        </a>
        . You can read the source code{" "}
        <a target="_blank" rel="noreferrer" href="https://github.com/Multiomics-Analytics-Group">
          on GitHub
        </a>
        .
      </p>
      <p>
        Nodes sizes are related to their{" "}
        <a target="_blank" rel="noreferrer" href="https://en.wikipedia.org/wiki/Betweenness_centrality">
          betweenness centrality
        </a>
        . More central nodes (ie. bigger nodes) are important crossing points in the network.
      </p>
    </Panel>
  );
};

export default DescriptionPanel;
