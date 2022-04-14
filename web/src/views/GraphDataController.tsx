//import Sigma from "sigma";
import { useSigma } from "react-sigma-v2";
import { FC, useEffect } from "react";
import { keyBy, omit } from "lodash";
import saveAsPNG from "./saveAsPNG";

import { Dataset, FiltersState } from "../types";

const GraphDataController: FC<{ dataset: Dataset; filters: FiltersState }> = ({ dataset, filters, children }) => {
  const sigma = useSigma();
  const graph = sigma.getGraph();
  const container = document.getElementById("sigma-container") as HTMLElement;
  /**const renderer = new Sigma(graph, container, {
    renderEdgeLabels: true,
  });**/
  
  /**
   * Feed graphology with the new dataset:
   */
  useEffect(() => {
    if (!graph || !dataset) return;

    const clusters = keyBy(dataset.clusters, "key");
    const tags = keyBy(dataset.tags, "key");

    dataset.nodes.forEach((node) => {
      graph.addNode(node.key, {
        ...node,
        ...omit(clusters[node.cluster], "key"),
        image: `${process.env.PUBLIC_URL}/images/${tags[node.tag].image}`,
      });}
    );
    dataset.edges.forEach((e) => { 
      const arr = Array.from(e);
      graph.addEdge(arr[0], arr[1], { size: arr[2] })});

    
    // Use degrees as node sizes:
    const scores = graph.nodes().map((node) => graph.getNodeAttribute(node, "score"));
    const minDegree = Math.min(...scores);
    const maxDegree = Math.max(...scores);
    const MIN_NODE_SIZE = 3;
    const MAX_NODE_SIZE = 30;
    graph.forEachNode((node) =>
      graph.setNodeAttribute(
        node,
        "size",
        ((graph.getNodeAttribute(node, "score") - minDegree) / (maxDegree - minDegree)) *
          (MAX_NODE_SIZE - MIN_NODE_SIZE) +
          MIN_NODE_SIZE,
      ),
    );

    graph.forEachNode((node) =>
      graph.setNodeAttribute(
        node,
        "hidden",
        false,)
    );



    
    return () => graph.clear();
  }, [graph, dataset]);


  

  /**
   * Apply filters to graphology:
   */
  useEffect(() => {
    const { clusters, tags } = filters;
    graph.forEachNode((node, { cluster, tag }) => 
      graph.setNodeAttribute(node, "hidden", !clusters[cluster] || !tags[tag]),
    );
    graph.forEachNode((node, {cluster, tag, hidden}) => {
      graph.forEachNeighbor(node, function(neighbor, attributes) {
        graph.setNodeAttribute(neighbor, "hidden", hidden || attributes['hidden']);
      });
    });
    
  }, [graph, filters]);


  // Bind save button:
const saveBtn = document.getElementById("save-as-png") as HTMLButtonElement;
/**saveBtn.addEventListener("click", () => {
  const layers = ["edges", "nodes", "edgeLabels", "labels"].filter(
    (id) => !!(document.getElementById(`layer-${id}`) as HTMLInputElement).checked,
  );

  saveAsPNG(renderer, layers);
});**/



  return <>{children}</>;
};

export default GraphDataController;
