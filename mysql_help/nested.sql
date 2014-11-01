CREATE Name Column WITH "ROW" TYPE AND "TABLE" TYPE
CREATE OR REPLACE TYPE "ROW_NAME"
  AS OBJECT
  (FIRST VARCHAR2(60),
   LAST VARCHAR2(60),
   MIDDLE VARCHAR2(30),
   PREFIX VARCHAR2(10),
   SUFFIX VARCHAR2(10),
   CREATE_DATE DATE,
   CHANGE_DATE DATE,
   CREATE_USER VARCHAR2(30),
   CHANGE_USER VARCHAR2(30))
/
CREATE OR REPLACE TYPE "TAB_NAME"
  AS TABLE OF ROW_NAME
/

--> Create Person Column with "ROW" Type and "TABLE" Type
CREATE OR REPLACE TYPE "ROW_PERSON"
  AS OBJECT
  (SSN VARCHAR2(9),
   GENDER VARCHAR2(1),
   BDATE DATE,
   ETHNICITY VARCHAR2(2),
   MARITAL_STATUS VARCHAR2(2),
   HAIR VARCHAR2(10),
   EYES VARCHAR2(10),
   HEIGHT VARCHAR2(10),
   CREATE_DATE DATE,
   CHANGE_DATE DATE,
   CREATE_USER VARCHAR2(30),
   CHANGE_USER VARCHAR2(30))
/
CREATE OR REPLACE TYPE "TAB_PERSON"
  AS TABLE OF ROW_PERSON
/

--> Create Employment Column with "ROW" Type and "TABLE" Type
CREATE OR REPLACE TYPE "ROW_EMPLOYMENT"
  AS OBJECT
  (JOB_TITLE VARCHAR2(60),
   START_DATE DATE,
   TERMINATED VARCHAR2(1),
   END_DATE DATE,
   SALARY NUMBER(8,4),
   CREATE_DATE DATE,
   CHANGE_DATE DATE,
   CREATE_USER VARCHAR2(30),
   CHANGE_USER VARCHAR2(30)
   )
/
CREATE OR REPLACE TYPE "TAB_EMPLOYMENT"
  AS TABLE OF ROW_EMPLOYMENT
/

--> Create table with UID and TAB_NAME & TAB_PERSON tables as columns
CREATE TABLE "ENTITY"
   (    "UIDN"          NUMBER(10) NOT NULL,
        "NAME"          "TAB_NAME" DEFAULT "TAB_NAME"(),
        "PERSON"        "TAB_PERSON" DEFAULT "TAB_PERSON"(),
        "EMPLOYMENT"    "TAB_EMPLOYMENT" DEFAULT "TAB_EMPLOYMENT"(),
        "DEAD"          VARCHAR2(1) DEFAULT 'N',
        "DEAD_DATE"     DATE,
        "CREATE_DATE"   DATE DEFAULT SYSDATE NOT NULL ENABLE,
        "CHANGE_DATE"   DATE,
        "CREATE_USER"   VARCHAR2(30) DEFAULT USER,
        "CHANGE_USER"   VARCHAR2(30),
         PRIMARY KEY ("UIDN"))
 NESTED TABLE "NAME" STORE AS "ENTITY_NAME"
 NESTED TABLE "PERSON" STORE AS "ENTITY_PERSON"
 NESTED TABLE "EMPLOYMENT" STORE AS "ENTITY_EMPLOYMENT"
/
COMMENT ON TABLE "ENTITY" IS 'Entity Base Table with repeating data as Nested Tables';
COMMENT ON COLUMN "ENTITY"."UIDN" IS 'Unique Identifaction Number/KEY';
COMMENT ON COLUMN "ENTITY"."NAME" IS 'Nested Table Storing Name Data for an entity';
COMMENT ON COLUMN "ENTITY"."PERSON" IS 'Nested Table Storing Person Data for an entity';
COMMENT ON COLUMN "ENTITY"."EMPLOYMENT" IS 'Nested Table Storing Employment related data for an entity';
COMMENT ON COLUMN "ENTITY"."DEAD" IS 'Y if entity is dead else N';
COMMENT ON COLUMN "ENTITY"."DEAD_DATE" IS 'Date of Death';
 
--> Now lets insert some records
INSERT INTO entity VALUES
(1,TAB_NAME(ROW_NAME('Orion','Pax',NULL,NULL,NULL,SYSDATE,NULL,USER,NULL),
            ROW_NAME('Optimus','Prime',NULL,NULL,NULL,SYSDATE,NULL,USER,NULL)),
 TAB_PERSON(ROW_PERSON('000000000','M','01-JAN-1985','AB',NULL,NULL,NULL,'10 meters',SYSDATE,NULL,USER,NULL)),
 TAB_EMPLOYMENT(ROW_EMPLOYMENT('Autobot Leader','01-JAN-1985',NULL,NULL,1000,SYSDATE,NULL,USER,NULL)),
 'N',NULL,SYSDATE,NULL,USER,NULL);
INSERT INTO entity VALUES
(2,TAB_NAME(ROW_NAME('Megatron',NULL,NULL,NULL,NULL,SYSDATE,NULL,USER,NULL)),
 TAB_PERSON(ROW_PERSON('000000001','M','01-JAN-1985','DC',NULL,NULL,NULL,'10 meters',SYSDATE,NULL,USER,NULL)),
 TAB_EMPLOYMENT(ROW_EMPLOYMENT('Decepticon Leader','01-JAN-1985',NULL,NULL,995,SYSDATE,NULL,USER,NULL)),
 'N',NULL,SYSDATE,NULL,USER,NULL);
COMMIT;
 
--> View what we have inserted by flattening out the column collections
SELECT src.uidn,
       n.*,
       p.*,
       e.*,
       dead,
       dead_date,
       src.create_date,
       src.change_date,
       src.create_user,
       src.change_user
FROM ENTITY src,
     TABLE(src.name) n,
     TABLE(src.person) p,
     TABLE(src.employment) e;
 
--> Insert specific records into column of collections
 
/* Megatron's name changed during Transformers the Movie so lets insert
   a new record into the collection for Names.
   The first thing that needs to be done is to get the base record from
   the entity and the column with the collection that we want to update
   table we do so with this statement
     INSERT INTO
     TABLE(select entity.name
           from entity
           where entity.uidn = 2)
   The next part of the insert involves adding a record to the collection
   column
     ('Galvatron',NULL,NULL,NULL,
       NULL,SYSDATE,NULL,USER,NULL)
*/
 
INSERT INTO
TABLE(SELECT entity.name
      FROM entity
      WHERE entity.uidn = 2)
VALUES
('Galvatron',NULL,NULL,NULL,
 NULL,SYSDATE,NULL,USER,NULL);
COMMIT;
 
--> Update specific records inside of column collections
/* We'll use a similiar method to update Megatrons Galvatron's
   NAME record to have a last name
*/
 
UPDATE TABLE(SELECT entity.name
             FROM entity
             WHERE entity.uidn = 2)
SET LAST = 'Decepticon',
    change_date = SYSDATE,
    change_user = USER
WHERE FIRST = 'Galvatron';
COMMIT;
 
--> Delete specific records inside of column collections
 
/* Eventually Galvatron's name returns to being Megatron,
   so lets delete out the Galvatron record from the NAME
   collection
*/

DELETE TABLE(SELECT entity.name
             FROM entity
             WHERE entity.uidn = 2)
WHERE FIRST = 'Galvatron';
COMMIT;