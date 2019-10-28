import java.util.*;


public class Player {
    private int PLAYER;
    private static final int RED_SIDE = 11;
    private static final int WHITE_SIDE = 20;
    private static int MAX_DEPTH = 9;
    private GameState bestState = new GameState();
    private HashMap<GameState, Integer> visitedStates = new HashMap<GameState, Integer>();

    public GameState play(final GameState gameState, final Deadline deadline) {
        Vector<GameState> nextStates = new Vector<GameState>();
        gameState.findPossibleMoves(nextStates);
        int size = nextStates.size();
        if (size == 0) {
            // Must play "pass" move if there are no other moves possible.
            return new GameState(gameState, new Move());
        }
        // Find out the current player
        PLAYER = gameState.getNextPlayer();
       
        for (int depth = 1; depth <= MAX_DEPTH; depth++) {
            alphaBeta(gameState, depth, Integer.MIN_VALUE, Integer.MAX_VALUE, PLAYER);
        }
        return bestState;
    }
    
    public int alphaBeta(GameState gameState, int depth, int alpha, int beta, int player) {
        
        if (gameState.isEOG() || depth == 0) {
            return getHeuristic(gameState);
        }
        Vector<GameState> nextStates = new Vector<>();
        GameState localBest = new GameState();
        gameState.findPossibleMoves(nextStates);
        int value;
        
        
        
        if (player == Constants.CELL_RED) {
            int maxValue = Integer.MIN_VALUE;
            value = Integer.MIN_VALUE;
            // Move ordering: sort the next moves by heuristic scores
            nextStates.sort(Comparator.comparingInt(this::getHeuristic));
            Collections.reverse(nextStates);
            for (int child = 0; child < nextStates.size(); child++) {
                GameState state = nextStates.get(child);
                // Avoid visiting repeated states
                if(visitedStates.containsKey(state)) {
                	value = visitedStates.get(state);  	
                }
                else {
                value = Math.max(value,
                        alphaBeta(
                                state, depth - 1, alpha, beta,Constants.CELL_WHITE));
                }
                alpha = Math.max(alpha, value);
                if (value > maxValue) {
                    maxValue = value;
                    localBest = state;
                }
                if (beta <= alpha)
                    break;      
            }
            bestState = localBest;
        }
        else {
            int minValue = Integer.MAX_VALUE;
            value = Integer.MAX_VALUE;
            // Move ordering: sort the next moves by heuristic scores
            nextStates.sort(Comparator.comparingInt(this::getHeuristic));
            for (int child = 0; child < nextStates.size(); child++) {
                GameState state = nextStates.get(child);
                // Avoid visiting repeated states
                if(visitedStates.containsKey(state)) {
                	value = visitedStates.get(state);  	
                }
                else {
                value = Math.min(value,
                        alphaBeta(
                                state, depth - 1, alpha, beta, Constants.CELL_RED));
                }
                beta = Math.min(beta, value);
                if (value < minValue) {
                    minValue = value;
                    localBest = state;
                }
                if (beta <= alpha)
                    break;     
            }
            bestState = localBest;
        }
       
        if(!visitedStates.containsKey(gameState)) visitedStates.put(gameState, value);
        
        return value;
    }
    
    
    
    public int getHeuristic(GameState gameState) {
        if (gameState.isEOG()) {
            if (gameState.isRedWin()) return Integer.MAX_VALUE;
            else if (gameState.isWhiteWin()) return Integer.MIN_VALUE;
            else if (gameState.isDraw()) return 0;
        }
        int whiteScore = 0, redScore = 0;
        for (int i = 0; i < GameState.NUMBER_OF_SQUARES; i++) {
            int piece = gameState.get(i);
            
            if ((piece & Constants.CELL_RED) != 0) {
                // King pieces worth higher values
                if ((piece & Constants.CELL_KING) != 0) {
                    redScore += 21;
                }else {
                    redScore += 1;
                }
                // A red piece reaching opponent's side
                if (i >= WHITE_SIDE) redScore += 7;
            }
            else if ((piece & Constants.CELL_WHITE) != 0) {
                // King pieces worth higher values
                if ((piece & Constants.CELL_KING) != 0) {
                    whiteScore += 21;
                } else {
                    whiteScore += 1;
                }
                // A red piece reaching opponent's side
                if (i <= RED_SIDE) whiteScore += 7;
            }
        }
        return redScore - whiteScore;
    }
}