
CREATE TABLE dxItems(
	item_id VARCHAR(100) NOT NULL,
	PRIMARY KEY(item_id)
) ;--#split#--

CREATE TABLE dxScoring_rules(
	item_id VARCHAR(100) NOT NULL,
	response VARCHAR(100) NOT NULL,
	item_score INTEGER NOT NULL, 
	
	PRIMARY KEY (item_id, response),
	FOREIGN KEY (item_id) REFERENCES dxItems(item_id) ON UPDATE CASCADE ON DELETE CASCADE
) ;--#split#--


CREATE TABLE dxBooklets(
	booklet_id VARCHAR(100) NOT NULL,
	
	PRIMARY KEY (booklet_id)	
) ;--#split#-- 



CREATE TABLE dxBooklet_design(
	booklet_id VARCHAR(100) NOT NULL,
	item_id VARCHAR(100) NOT NULL,
	item_position INTEGER NOT NULL CHECK(item_position >= 1),
	
	PRIMARY KEY (booklet_id, item_id),
	UNIQUE		(booklet_id, item_position),
	
	FOREIGN KEY (item_id) REFERENCES dxItems(item_id) ON UPDATE CASCADE ON DELETE CASCADE,
	FOREIGN KEY (booklet_id) REFERENCES dxBooklets(booklet_id) ON UPDATE CASCADE ON DELETE CASCADE
) ;--#split#--


 
CREATE TABLE dxPersons(
	person_id VARCHAR(100) NOT NULL PRIMARY KEY
)  ;--#split#--



CREATE TABLE dxAdministrations(
	person_id VARCHAR(100) NOT NULL,
	booklet_id VARCHAR(100) NOT NULL,
	
	PRIMARY KEY (person_id,booklet_id),

	FOREIGN KEY (person_id) REFERENCES dxPersons(person_id) ON UPDATE CASCADE ON DELETE CASCADE,
	FOREIGN KEY (booklet_id) REFERENCES dxBooklets(booklet_id) ON UPDATE CASCADE ON DELETE CASCADE
) ;--#split#--

CREATE TABLE dxResponses(
	person_id VARCHAR(100) NOT NULL,
	booklet_id VARCHAR(100) NOT NULL,
	item_id VARCHAR(100) NOT NULL,
	response VARCHAR(100) NOT NULL,
	
	PRIMARY KEY (person_id, booklet_id, item_id),
	
	FOREIGN KEY (person_id,booklet_id) REFERENCES dxAdministrations(person_id,booklet_id) ON UPDATE CASCADE ON DELETE CASCADE DEFERRABLE INITIALLY DEFERRED,
	FOREIGN KEY (booklet_id, item_id) REFERENCES dxBooklet_design(booklet_id, item_id) ON UPDATE CASCADE ON DELETE CASCADE DEFERRABLE INITIALLY DEFERRED,
	FOREIGN KEY (item_id, response) REFERENCES dxScoring_rules(item_id, response) ON UPDATE CASCADE ON DELETE CASCADE DEFERRABLE INITIALLY DEFERRED
	-- foreign key constraints deferred for speed on some systems
) ;--#split#--

CREATE VIEW dxBooklet_stats AS
WITH 
A1 AS ( SELECT booklet_id, COUNT(*) AS n_persons FROM dxAdministrations GROUP BY booklet_id),
I1 AS ( SELECT booklet_id, COUNT(DISTINCT item_id) AS n_items FROM  dxBooklet_design  GROUP BY booklet_id)
SELECT B.booklet_id, COALESCE(n_items,0) AS n_items, COALESCE(n_persons,0) AS n_persons
	FROM dxBooklets AS B
		LEFT OUTER JOIN A1 USING(booklet_id)
			LEFT OUTER JOIN I1 USING(booklet_id);



