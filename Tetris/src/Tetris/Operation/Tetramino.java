package Tetris.Operation;

import java.awt.*;

/**
* Egy tetramino-t képvisel
 */
public abstract class Tetramino {

    private final Point position;
    private int[][] shape;
    private Color color;

    /**
     * Létrehoz egy új tetramino-t és beállítja a kezdő pozícióját
     */
    public Tetramino() {
        position = new Point(5,1);
        color = Color.WHITE;
    }

    /**
     * Megváltoztatja a tetramino alakzatát
     * @param newShape új alakzat
     */
    protected void changeShape(int[][] newShape){
        shape = newShape;
    }

    /**
     * Megadja a tetramino pozicióját
     * @return pozíció
     */
    public Point getPosition(){
        return position;
    }

    /**
     * Függőlegesen mozgatja a tetramino-t
     */
    public void dropTetramino(){
        (position.y)++;
    }

    /**
     * Mozgatja vízszintes irányba a tetramino-t
     * @param deltaX mozgatás mértéke
     */
    public void move(int deltaX){
        position.x+=deltaX;
    }

    /**
     * Megadja a tetramino alakzatát
     * @return alakzat
     */
    public int[][] getShape(){
        return shape;
    }

    /**
     * Beállítja a tetramino színét
     * @param c szín
     */
    public void setColor(Color c){
        color = c;
    }

    /**
     * Megadja a tetramino színét
     * @return tetramino színe
     */
    public Color getColor(){
        return color;
    }

    /**
     * Megad egy ugyanolyan tetramino-t
     * @return lemásolt tetramino
     */
    abstract public Tetramino copy();

}
