#ifndef SRC_DECOMPOSITION_ABSTRACT_STRATEGY_MODELS_HPP
#define SRC_DECOMPOSITION_ABSTRACT_STRATEGY_MODELS_HPP

using val_t = double;

struct ModelComputationalCosts {
  template <typename particles, typename DecompositionStrategy,
            typename DistributionStrategy>
  void calculate(particles &vd, DecompositionStrategy &dec,
                 DistributionStrategy &dist){};

  /*! \brief Calculate communication and migration costs
   *
   * \param ts how many timesteps have passed since last calculation, used to
   * approximate the cost
   */
  template <typename DecompositionStrategy, typename DistributionStrategy>
  void computeCommunicationAndMigrationCosts(DecompositionStrategy &dec,
                                             DistributionStrategy &dist,
                                             const size_t ts = 1){};
};

struct ModelDecompose {
  template <typename Decomposition>
  void applyModel(Decomposition &dec, size_t v){};

  template <typename DecompositionStrategy, typename DistributionStrategy>
  void finalize(DecompositionStrategy &dec, DistributionStrategy &dist){};
};

struct ModelDistribute {
  virtual val_t toll() = 0;
};

#endif // SRC_DECOMPOSITION_ABSTRACT_STRATEGY_MODELS_HPP
