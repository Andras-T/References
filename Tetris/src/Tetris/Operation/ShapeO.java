package Tetris.Operation;

/**
 * 'O' alakú tetramino
 */
public class ShapeO extends Tetramino{

    /**
     * Létrehoz egy 'O' alakú tetramino-t
     */
    public ShapeO(){
        super();
        int[][] shape = new int[][]{
                {0, 0, 0, 0},
                {0, 1, 1, 0},
                {0, 1, 1, 0},
                {0, 0, 0, 0}};
        changeShape(shape);
    }

    /**
     * Létrehoz egy másolatot
     * @return másolat
     */
    public ShapeO copy(){
        return new ShapeO();
    }

}
