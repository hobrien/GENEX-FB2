DROP TABLE IF EXISTS `GTEx_eQTLs`;
CREATE TABLE `GTEx_eQTLs` (
  `variant_id` varchar(255) NOT NULL,
  `gene_id` varchar(255) NOT NULL,
  `tss_distance` int(11) DEFAULT NULL,
  `minor_allele_samples` int(11) DEFAULT NULL,
  `minor_allele_count` int(11) DEFAULT NULL,
  `maf` float DEFAULT NULL,
  `pval_nominal` float DEFAULT NULL,
  `slope` float DEFAULT NULL,
  `slope_se` float DEFAULT NULL,
  `pval_nominal_threshold` float DEFAULT NULL,
  `min_pval_nominal` float DEFAULT NULL,
  `pval_beta` float DEFAULT NULL,
  `tissue` varchar(255) NOT NULL,
  PRIMARY KEY (`tissue`,`gene_id`,`variant_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;
