/*******************************************************************************
 * Copyright (c) 2013 Christian Wiwie.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Public License v3.0
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/gpl.html
 * 
 * Contributors:
 *     Christian Wiwie - initial API and implementation
 ******************************************************************************/
/**
 * 
 */
package de.clusteval.cluster.quality;

import java.io.File;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import de.clusteval.cluster.ICluster;
import de.clusteval.cluster.IClusterItem;
import de.clusteval.cluster.IClustering;
import de.clusteval.data.IDataConfig;
import de.clusteval.framework.repository.IRepository;
import de.clusteval.framework.repository.RegisterException;
import de.clusteval.utils.ClustEvalAlias;
import de.clusteval.utils.ClustEvalVersionRequirement;
import de.clusteval.utils.DynamicComponentVersion;


/**
 * @author Christian Wiwie
 */
@DynamicComponentVersion(version = "3-SNAPSHOT")
@ClustEvalAlias(alias = "F1-Score")
@ClustEvalVersionRequirement(target = IClusteringQualityMeasure.class, versionSpecification = "[2]")
public class TransClustFClusteringQualityMeasure extends ClusteringQualityMeasure {

	/**
	 * 
	 */
	private static final long serialVersionUID = -8568282992888457153L;

	/**
	 * @param repo
	 * @param changeDate
	 * @param absPath
	 * @throws RegisterException
	 */
	public TransClustFClusteringQualityMeasure(IRepository repo, long changeDate, File absPath,
			ClusteringQualityMeasureParameters parameters) throws RegisterException {
		super(repo, changeDate, absPath, parameters);
	}

	/**
	 * The copy constructor for this measure.
	 * 
	 * @param other
	 *            The object to clone.
	 * @throws RegisterException
	 */
	public TransClustFClusteringQualityMeasure(final TransClustFClusteringQualityMeasure other)
			throws RegisterException {
		super(other);
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see
	 * cluster.quality.ClusteringQualityMeasure#getQualityOfClustering(cluster
	 * .IClustering, data.goldstandard.GoldStandard)
	 */
	@SuppressWarnings("unused")
	@Override
	public ClusteringQualityMeasureValue getQualityOfClustering(final IClustering clustering, IClustering gsClustering,
			final IDataConfig dataConfig) {

		double fmeasure = 0;

		/*
		 * Ensure, that goldstandard contains only objects, that are also in the
		 * dataset. Otherwise precision will be calculated incorrectly, because
		 * it directly depends on the number of items in a cluster in the
		 * goldstandard.
		 */
		Set<IClusterItem> gsIClusterItems = new HashSet<IClusterItem>(gsClustering.getClusterItems());
		Set<IClusterItem> clusterItems = new HashSet<IClusterItem>(clustering.getClusterItems());
		gsIClusterItems.removeAll(clusterItems);
		for (IClusterItem onlyInGs : gsIClusterItems)
			gsClustering.removeClusterItem(onlyInGs);

		/*
		 * Ensure, that clustering contains only objects, that are also in the
		 * goldstandard.
		 */
		gsIClusterItems = new HashSet<IClusterItem>(gsClustering.getClusterItems());
		clusterItems.removeAll(gsIClusterItems);

		for (IClusterItem onlyInClustering : clusterItems)
			clustering.removeClusterItem(onlyInClustering);

		final float proteins = clustering.fuzzySize();

		for (ICluster gsCluster : gsClustering) {
			final float proteinsInReference = gsCluster.fuzzySize();
			// final double maxValue = findMax(clustering, gsCluster);
			final double maxValue = findMax2(clustering, gsCluster);
			fmeasure += (maxValue * proteinsInReference);
		}
		fmeasure /= proteins;

		return ClusteringQualityMeasureValue.getForDouble(fmeasure);
	}

	/**
	 * Find max.
	 * 
	 * @param proteinsInReference
	 *            the proteins in reference
	 * @param clustering
	 *            the clustering
	 * @param gsCluster
	 *            the gs cluster
	 * @return the double
	 */
	private static double findMax2(final IClustering clustering, final ICluster gsCluster) {
		double max = 0;
		double maxCommon = 0;

		for (ICluster cluster : clustering) {
			double common = 0;

			// performance reasons
			if (gsCluster.size() < cluster.size()) {
				common = calculateCommonProteins(gsCluster, cluster);
			} else {
				common = calculateCommonProteins(cluster, gsCluster);
			}

			final double tp = common;
			final double fp = cluster.fuzzySize() - common;
			final double fn = gsCluster.fuzzySize() - common;
			final double precision = (tp / (tp + fp));
			final double recall = (tp / (tp + fn));
			final double fmeasure = 2 * (recall * precision) / (recall + precision);

			if (fmeasure > max) {
				max = fmeasure;
				maxCommon = common;
			}
		}
		if (maxCommon == 0 && gsCluster.size() == 1) {
			return 1;
		}
		return max;
	}

	/**
	 * Calculate common proteins.
	 * 
	 * @param c1
	 *            the c1
	 * @param c2
	 *            the c2
	 * @return the float
	 */
	private static float calculateCommonProteins(final ICluster c1, final ICluster c2) {
		float common = 0;

		Map<IClusterItem, Float> items = c1.getFuzzyItems();

		for (IClusterItem item : items.keySet()) {
			if (c2.contains(item)) {
				common += items.get(item);
			}
		}
		return common;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see cluster.quality.ClusteringQualityMeasure#getMinimum()
	 */
	@Override
	public double getMinimum() {
		return 0.0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see cluster.quality.ClusteringQualityMeasure#getMaximum()
	 */
	@Override
	public double getMaximum() {
		return 1.0;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see cluster.quality.ClusteringQualityMeasure#requiresGoldstandard()
	 */
	@Override
	public boolean requiresGoldstandard() {
		return true;
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see cluster.quality.ClusteringQualityMeasure#isBetterThanHelper(cluster.
	 * quality .ClusteringQualityMeasureValue,
	 * cluster.quality.ClusteringQualityMeasureValue)
	 */
	@Override
	protected boolean isBetterThanHelper(IClusteringQualityMeasureValue quality1,
			IClusteringQualityMeasureValue quality2) {
		return quality1.getValue() > quality2.getValue();
	}

	/*
	 * (non-Javadoc)
	 * 
	 * @see de.clusteval.cluster.quality.ClusteringQualityMeasure#
	 * supportsFuzzyClusterings()
	 */
	@Override
	public boolean supportsFuzzyClusterings() {
		return false;
	}
}

