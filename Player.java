import java.util.*;
import java.util.stream.Collectors;

class Player {
    // Each kind of bird have several models
    private int N;      // number of states in HMM
    private int M;      // number of emissions in HMM
    private ArrayList<HMM>[] speciesHMMs;   // Store some models for each specie
    private int[] lGuess;       // the guess I made
    private int timeStep;       // the time elapsed in each round
    private int currentRound;   // record the current round
    private int correctGuess;   // number of correct guesses in each round
    private int wrongGuess;     // number of wrong guesses in each round
    private int totalCorrectGuess;  // total number of correct guesses
    private int totalWrongGuess;    // total number of wrong guesses
    private int shootTime;          // the time allowed to start shooting
    private int hits;               // record how many birds we hit
    private int triedShoots;        // record the number of trials
    private double latestCorrectness; // record the accuracy in the last round

    class GuessResult {
        private final int species;
        private final int model;
        private final double prob;

        public GuessResult(int species, int model, double prob) {
            this.species = species;
            this.model = model;
            this.prob = prob;
        }
    }
    public Player() {
        shootTime = 20;
        N = Constants.COUNT_SPECIES;
        M = Constants.COUNT_MOVE;
        speciesHMMs = new ArrayList[Constants.COUNT_SPECIES];
        for (int i = 0; i < N; i++) {
            speciesHMMs[i] = new ArrayList<>();
            speciesHMMs[i].add(new HMM(N, M));
        }
    }

    /**
     * Shoot!
     *
     * This is the function where you start your work.
     *
     * You will receive a variable pState, which contains information about all
     * birds, both dead and alive. Each bird contains all past moves.
     *
     * The state also contains the scores for all players and the number of
     * time steps elapsed since the last time this function was called.
     *
     * @param pState the GameState object with observations etc
     * @param pDue time before which we must have returned
     * @return the prediction of a bird we want to shoot at, or cDontShoot to pass
     */
    public Action shoot(GameState pState, Deadline pDue) {
        /*
         * Here you should write your clever algorithms to get the best action.
         * This skeleton never shoots.
         */
        timeStep++;
        if (pState.getRound() > currentRound) {
            // Reset the statics in each new round
            correctGuess = 0;
            wrongGuess = 0;
            timeStep = 0;
            currentRound = pState.getRound();
        }
        if (timeStep < 100 - shootTime) {
            // don' t shoot before we get long enough observations
            return cDontShoot;
        }
        if (currentRound > 0) {
            // don't shoot in the 1st round, not prepared
            return decideShootAction(pState);
        }
        return cDontShoot;
    }

    public Action decideShootAction(GameState pState) {
        double globalMaxProb = 0.;
        int chosenBird = -1;
        int globalBestAction = -1;
        for (int i = 0; i < pState.getNumBirds(); i++) {
            Bird bird = pState.getBird(i);
            // Don't shoot dead birds
            if (bird.isDead()) continue;
            Integer[] obs = getObs(bird);
            int label = conservativeGuess(obs).species;
            // Don't shoot black storks
            if (label == Constants.SPECIES_BLACK_STORK) {
                continue;
            }
            // Train a new model for the current bird
            HMM model = new HMM(N - 1, M);
            model.BWAlgorithm(obs);
            /* Apply forward algorithm to obtain the probability of being at each hidden state
             at time T */
            double[][] alpha = model.getAlpha();
            double[] stateDistribution = Matrices.getMatCol(alpha, model.getT() - 1);
            Matrices.normalize(stateDistribution);
            double[] nextEmission = model.nextEmission(stateDistribution);
            Matrices.normalize(nextEmission);

            double max = 0.;
            int mostLikelyMove = 0;
            for (int j = 0; j < nextEmission.length; j++) {
                if (nextEmission[j] > max) {
                    max = nextEmission[j];
                    mostLikelyMove = j;
                }
            }
            // Update the overall optimal choice
            if (max > globalMaxProb) {
                globalMaxProb = max;
                globalBestAction = mostLikelyMove;
            }
        }
        // Set a threshold for shooting to guarantee the hit rates
        if (globalMaxProb > 0.98) {
            this.triedShoots++;
            System.err.println("Shoot with prob: " + globalMaxProb);
            return new Action(chosenBird, globalBestAction);
        }else {
            return cDontShoot;
        }
    }

    /**
     * Guess the species!
     * This function will be called at the end of each round, to give you
     * a chance to identify the species of the birds for extra points.
     *
     * Fill the vector with guesses for the all birds.
     * Use SPECIES_UNKNOWN to avoid guessing.
     *
     * @param pState the GameState object with observations etc
     * @param pDue time before which we must have returned
     * @return a vector with guesses for all the birds
     */
    public int[] guess(GameState pState, Deadline pDue) {
        /*
         * Here you should write your clever algorithms to guess the species of
         * each bird. This skeleton makes no guesses, better safe than sorry!
         */
        printInfo(pState);
        Random random = new Random();
        this.lGuess = new int[pState.getNumBirds()];
        if (pState.getRound() == 0) {
            for (int i = 0; i < pState.getNumBirds(); ++i)
                lGuess[i] = 0;
        } else {
            for (int i = 0 ; i < pState.getNumBirds(); i++) {
                // Extract observations sequence for each bird
                Bird bird = pState.getBird(i);
                Integer[] obs = getObs(bird);
                if (obs.length == 0) {
                    lGuess[i] = random.nextInt(N - 1);
                } else {
                    GuessResult myGuess = conservativeGuess(obs);
                    lGuess[i] = myGuess.species;
                }
            }
        }
        System.err.println("My guess");
        for (int i = 0; i < lGuess.length; i++) {
            System.err.print(lGuess[i] + " ");
        }
        System.err.println();
        return lGuess;
    }

    public GuessResult conservativeGuess(Integer[] obs) {
        // Record how likely the bird is of that species
        double[] speciesProbs = new double[N];
        // The sampling probability
        double chooseProb = 0.98;
        // For each species, find the model thatfits the observation most
        int[] bestModels = new int[N];
        // Traverse all the models
        for (int i = 0; i < N; i++) {
            double max = 0.;
            int model = 0;
            for (int j = 0; j < speciesHMMs[i].size(); j++){
                // Randomly pick some models
                if (Math.random() > chooseProb && speciesHMMs[i].size() > 100) {
                    continue;
                }
                // Estimate how likely the model would yield such sequences
                double prob = speciesHMMs[i].get(j).forwardAlgorithm(obs);
                if (prob > max) {
                    max = prob;
                    model = j;
                }
            }
            speciesProbs[i] = max;
            bestModels[i] = model;
        }
        Matrices.normalize(speciesProbs);
        // Decide the species of the bird according to the likelihood
        int bestModel = 0;
        int bestSpecies = 0;
        double bestProb = 0.;
        for (int i = 0; i < N; i++) {
            if (speciesProbs[i] > bestProb) {
                bestProb = speciesProbs[i];
                bestSpecies = i;
                bestModel = bestModels[i];
            }
        }
        // P < 0.5 means it could be a black stork, because right guesses always yield 0.99 or more
        if (bestProb < 0.5) {
            return new GuessResult(Constants.SPECIES_BLACK_STORK, -1, 0.);
        }
        return new GuessResult(bestSpecies, bestModel, bestProb);
    }

    /*
        Get the observation sequences for given bird
     */
    public Integer[] getObs(Bird bird) {
        ArrayList<Integer> list = new ArrayList<>();
        for (int i = 0; i < bird.getSeqLength(); i++) {
            if (bird.wasDead(i)) {
                break;
            }
            list.add(bird.getObservation(i));
        }
        Integer[] obs = list.toArray(new Integer[list.size()]);
        return obs;
    }

    public void printInfo(GameState pState) {
        System.err.println("************ Round " + pState.getRound() +" ************");
        System.err.println("Number of birds: " + pState.getNumBirds());
    }

    /**
     * If you hit the bird you were trying to shoot, you will be notified
     * through this function.
     *
     * @param pState the GameState object with observations etc
     * @param pBird the bird you hit
     * @param pDue time before which we must have returned
     */
    public void hit(GameState pState, int pBird, Deadline pDue) {
        hits++;
        System.err.println("HIT BIRD!!!");
    }

    /**
     * If you made any guesses, you will find out the true species of those
     * birds through this function.
     *
     * @param pState the GameState object with observations etc
     * @param pSpecies the vector with species
     * @param pDue time before which we must have returned
     */
    public void reveal(GameState pState, int[] pSpecies, Deadline pDue) {
        // Store the i-th observation sequences
        ArrayList<Integer[]> roundObs = new ArrayList<>();
        for (int i = 0; i < pState.getNumBirds(); i++) {
            roundObs.add(getObs(pState.getBird(i)));
        }
        // Collect all the observations sequences for each species respectively
        HashMap<Integer, ArrayList<Integer[]>> sortedObs = new HashMap<>();
        for (int i = 0; i < Constants.COUNT_SPECIES; i++) {
            sortedObs.put(i, new ArrayList<>());
        }
        for (int i = 0; i < pSpecies.length; i++) {
            sortedObs.get(pSpecies[i]).add(roundObs.get(i));
        }
        // Train one model with multiple observation sequences
        for (int i = 0; i < Constants.COUNT_SPECIES; i++) {
            ArrayList<Integer[]> specieObs = sortedObs.get(i);
            if (specieObs.size() == 0) continue;
            // Repeat several times to avoid local optimum
            for (int k = 0; k < 8; k++) {
                // Avoid overfitting
                HMM specialModel = new HMM(1, M, 20);
                specialModel.trainMultiSeq(specieObs);
                speciesHMMs[i].add(specialModel);
            }
        }
        // Simple statics
        for (int i= 0; i < pSpecies.length; i++) {
            if (pSpecies[i] == lGuess[i]) {
                correctGuess++;
                totalCorrectGuess++;
            }
            else {
                wrongGuess++;
                totalWrongGuess++;
            }
        }
        // Display the information
        this.latestCorrectness = (double) correctGuess / (correctGuess + wrongGuess);
        System.err.println("In this round: correct guesses:" + correctGuess + ", wrong guesses:" +
                wrongGuess + " Accuracy: " + latestCorrectness);
        System.err.println("Totally: correct guesses:" + totalCorrectGuess + ", wrong guesses:" +
                totalWrongGuess + " Accuracy: " + (double) totalCorrectGuess / (totalCorrectGuess + totalWrongGuess));
        System.err.println("Hits/Trials:" + hits + " " + this.triedShoots);
        System.err.println("Review answers");
        for (int i = 0; i < pSpecies.length; i++) {
            System.err.print(pSpecies[i] + " ");
        }
        System.err.println();
    }
    public static final Action cDontShoot = new Action(-1, -1);
}
