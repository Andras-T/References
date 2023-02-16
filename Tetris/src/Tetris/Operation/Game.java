package Tetris.Operation;

import Tetris.Graphics.Board;
import javax.swing.*;
import java.awt.*;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Random;

/**
 * Tetris játék
 */
public class Game implements KeyListener{

    //tárolók
    private final ArrayList<Tetramino> nextPieces = new ArrayList<>();
    private final ArrayList<Color> nextColor = new ArrayList<>();

    // Képfrissítés
    private static final int FPS = 30;
    private static final int delay = FPS/1000;
    private Timer looper;

    // Sebesség
    private int normal = 450;
    private static final int fast = 100;
    private static final int drop = 10;
    private int delayTimeForMovement = normal;
    private static long beginTime;

    // Pálya
    private static final int WIDTH = 14, HEIGHT = 24;
    private Color[][] well;
    private Tetramino currentPiece;
    private int score;
    private boolean paused;
    private final Board board;
    protected final GameWindow gameWindow;

    /**
     * Létrehoz egy új játékot
     * @param window az ablak amelyen a játék fut
     */
    public Game(GameWindow window){
        gameWindow = window;
        board = window.getBoard();
        board.addKeyListener(this);
        initialize();
        loop();
    }

    /**
     * Inicializálja az adatokat
     */
    protected void initialize(){
        paused = false;

        // Tetraminok hozzáadása
        nextPieces.add(new ShapeI());
        nextPieces.add(new ShapeZ());
        nextPieces.add(new ShapeS());
        nextPieces.add(new ShapeJ());
        nextPieces.add(new ShapeL());
        nextPieces.add(new ShapeT());
        nextPieces.add(new ShapeO());
        // Színek hozzáadása
        nextColor.add(Color.cyan);
        nextColor.add(Color.blue);
        nextColor.add(Color.red);
        nextColor.add(Color.pink);
        nextColor.add(Color.yellow);
        nextColor.add(Color.orange);
        nextColor.add(Color.green);
        nextColor.add(Color.white);
        nextColor.add(Color.MAGENTA);

        //tábla feltöltése
        well = new Color[WIDTH][HEIGHT];
        for (int i = 0; i < WIDTH; i++) {
            for (int j = 0; j < HEIGHT; j++) {
                if (i == 0 || i == 1 || i == 12 || i == 13 || j == 22 || j == 23) {
                    well[i][j] = Color.GRAY;
                } else {
                    well[i][j] = Color.BLACK;
                }
            }
        }
        newPiece();
    }

    /**
     * Új tetramino-t sorsol
     */
    public void newPiece() {
        Collections.shuffle(nextPieces, new Random());
        Collections.shuffle(nextColor,  new Random());
        currentPiece = (nextPieces.get(2)).copy();
        currentPiece.setColor(nextColor.get(1));
    }

    /**
     * Frissíti a játékot
     */
    private void loop(){
        looper = new Timer(delay, event -> {
            paintWell();

            // Függőleges mozgás
            if (System.currentTimeMillis() - beginTime >= delayTimeForMovement) {
                if (!paused) {
                    //Gyorsítás
                    if (normal > 250)
                        normal--;
                    drop();
                }
                beginTime = System.currentTimeMillis();
            }
        });
        looper.start();
    }

    /**
     * Oldal irányú mozgás
     * @param deltaX mozgás iránya és távolsága
     */
    public void move(int deltaX){
        if (!collides(deltaX,0))
            currentPiece.move(deltaX);
    }

    /**
     * Függőleges irányú egységnyi mozgás
     */
    public void drop(){
        if (!collides(0,1)) {
            currentPiece.dropTetramino();
        }
        else{
            if (currentPiece.getPosition().y <= 1) {
                gameOver();
            }
            fixToWell();
            newPiece();
        }
    }

    /**
     * Teszteli, hogy ütközik-e a tetramino egy adott irányú mozgás esetén
     * @param deltaX Vízszintes irányú mozgás
     * @param deltaY Függőleges irányú mozgás
     * @return ütközik-e
     */
    public boolean collides(int deltaX, int deltaY){
        int[][] shape = currentPiece.getShape();
        int length = shape.length;
        Point position = currentPiece.getPosition();
        for (int y = 0; y < length; y++){
            for (int x = 0; x < length; x++){
                if(shape[x][y] != 0) {
                    if (well[position.x + x +deltaX][position.y + y + deltaY] != Color.BLACK){
                        return true;
                    }
                }
            }
        }
        return false;
    }

    /**
     * Forgatja az aktuális tetramino-t
     * @param current Aktuális tetramino
     */
    public void rotate(Tetramino current){
        Point position = current.getPosition();
        int[][] shape = current.getShape();
        // Alakzat forgatása
        try {
            int[][] rotatedShape = rotateCW(shape);

        //alakzat alkalmazása a tetramino-n
        for (int y = 0; y < 4; y++) {
            for (int x = 0; x < 4; x++) {
                if (rotatedShape[x][y] != 0) {
                    if(well[position.x + x][position.y + y] != Color.BLACK)
                        return;
                }
            }
        }
        current.changeShape(rotatedShape);
        }catch (IndexOutOfBoundsException | NullPointerException e){
            e.printStackTrace();
        }
    }

    /**
     * Óramutató járásával megegyező irányba forgat egy NxN-es alakzatot
     * @param shape elforgatni kívánt alakzat
     * @return  elforgatott alakzat
     * @throws IndexOutOfBoundsException Nem megfelelő méretű alakzat esetén
     * @throws NullPointerException Üres alakzat esetén
     */
     public int[][] rotateCW(int[][] shape) throws IndexOutOfBoundsException, NullPointerException{
        final int length = shape[0].length;
        int[][] rotatedShape = new int[length][length];
        for (int row = 0; row < length; row++) {
            for (int column = 0; column < length; column++) {
                rotatedShape[column][length-1-row] = shape[row][column];
            }
        }
        return rotatedShape;
    }

    /**
     * Rögzít egy tetramino-t a táblához
     */
    public void fixToWell() {
        int[][] shape = currentPiece.getShape();
        final int length = shape.length;
        Color tetraminoColor = currentPiece.getColor();
        Point position = currentPiece.getPosition();

        for (int row = 0; row < length; row++){
            for (int column = 0; column < length; column++){
                if (shape[column][row] != 0) {
                    well[position.x + column][position.y + row] = tetraminoColor;
                }
            }
        }
        searchFullRow();
    }

    /**
     * Teli sort keres a pályán és ha talál törli és a felette lévő sorokat lejjebb helyezi
     */
    public void searchFullRow(){
        for (int row = HEIGHT-3; row > 1; row--) {
            int db = 0;
            for (int column = 2;column < WIDTH-2; column++) {
                if (well[column][row] != Color.BLACK)
                    db++;
                if (db == 10){
                    score+=10;
                    deleteRow(row);
                    searchFullRow();
                }
            }
        }
    }

    /**
     * Töröl egy sort a pályáról és a felette lévő sorokat lejjebb helyezi
     * @param delete Az adott sor sorszáma
     */
    public void deleteRow(int delete){
        for (int row = delete; row > 1 ; row--){
            for (int column = 2; column < WIDTH-2; column++)
            {
                well[column][row] = well[column][row-1];
            }
        }
    }

    /**
     * Kirajzolja a táblát a pályára
     */
    public void paintWell()
    {
        Color[][] newWell = new Color[WIDTH][HEIGHT];
        for (int i = 0; i < 14; i++) {
            newWell[i] = well[i].clone();
        }

        int[][] shape = currentPiece.getShape();
        final int length = shape.length;
        Color tetraminoColor = currentPiece.getColor();
        Point position = currentPiece.getPosition();

        for (int y = 0; y < length; y++){
            for (int x = 0; x < length; x++){
                if (shape[x][y] != 0) {
                    newWell[position.x + x][position.y + y] = tetraminoColor;
                }
            }
        }
        board.paintWell(newWell);
    }

    /**
     * Játék vége
     */
    public void gameOver(){
        looper.stop();
        try {
            save(score);
        } catch (IOException e) {
            e.printStackTrace();
        }
        gameWindow.gameOver(score);
    }

    /**
     * Eredmény mentése, ha a legjobb 5 közé sorolható
     * @throws IOException Nem sikerült megnyitni az eredményekhez tartozó fájlt
     */
    public void save(int score)throws IOException{
        //beolvasás
        FileReader fr = new FileReader("./res/HighScores.txt");
        BufferedReader br = new BufferedReader(fr);
        Integer[] scores = new Integer[6];
        for (int i = 0; i < 5; i++){
            String line = br.readLine();
            if (line == null)
            {
                scores[i] = 0;
            }
            else
            scores[i] = Integer.parseInt(line);
        }
        br.close();
        scores[5] = score;
        Arrays.sort(scores, Collections.reverseOrder());

        //kiiratás
        FileWriter fw = new FileWriter("./res/HighScores.txt");
        PrintWriter pw = new PrintWriter(fw);
        for (int i = 0; i < 5; i++){
            pw.println(scores[i]);
        }
        pw.close();
    }

    @Override
    public void keyTyped(KeyEvent e) {}

    @Override
    public void keyPressed(KeyEvent e) {
        switch (e.getKeyCode()) {

            case KeyEvent.VK_UP:
                if(!paused)
               rotate(currentPiece);
                break;
            case KeyEvent.VK_DOWN:
                if(!paused)
                delayTimeForMovement = fast;
                break;
            case KeyEvent.VK_LEFT:
                if(!paused)
                move(-1);
                break;
            case KeyEvent.VK_RIGHT:
                if(!paused)
                move(+1);
                break;
            case KeyEvent.VK_SPACE:
                if(!paused)
                delayTimeForMovement = drop;
                break;
            case KeyEvent.VK_ESCAPE:
                System.out.println(e.getKeyChar());
                paused= !paused;
                board.pausePanel(paused);
                break;
        }
    }

    @Override
    public void keyReleased(KeyEvent e) {
        if (e.getKeyCode() == KeyEvent.VK_SPACE)
            delayTimeForMovement = normal;
        else if (e.getKeyCode() == KeyEvent.VK_DOWN)
            delayTimeForMovement = normal;
    }
}
