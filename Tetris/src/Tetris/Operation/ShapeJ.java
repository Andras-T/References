package Tetris.Operation;

/**
 * 'J' alakú tetramino
 */
public class ShapeJ extends Tetramino{

    /**
     * Létrehoz egy 'J' alakú tetramino-t
     */
    public ShapeJ(){
        super();
        int[][] shape = new int[][]{
                {0, 0, 1, 0},
                {0, 0, 1, 0},
                {0, 1, 1, 0},
                {0, 0, 0, 0}};
        changeShape(shape);
    }

    /**
     * Létrehoz egy másolatot
     * @return másolat
     */
    public ShapeJ copy(){
        return new ShapeJ();
    }

}
